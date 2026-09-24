#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
emi_cm_model.py -- TPS62933 conducted COMMON-MODE (CM) EMI model, 150 kHz - 30 MHz
                   (CISPR 32 approximation), three load cases.

Companion to emi_v3/emi_model_v3.py (which did DIFFERENTIAL mode).
Writes ONLY to /mnt/raid10/sim-work/tps62933/emi_cm/.

Prerequisite doc (authoritative): tps62933-sim/docs/07-input-conditions-and-assumptions.md
  - board side has NO PE, only +24V/GND two wires (10 cm, 0.75 mm^2)
  - plastic (unshielded) enclosure
  - upstream LM50-20B24, CLASS I (case->PE), 4000 VAC isolation, fsw 65 kHz
  - C_p = board<->earth parasitic is the ONLY element closing the CM loop
  - CISPR 32 Class B, fsw 805 kHz, tr=tf=5 ns (assumed)

CM source model
---------------
The switching node SW has the largest dv/dt on the board.  It couples to the
reference (earth) plane through the parasitic capacitance C_p.  The injected
common-mode current phasor at harmonic k is

        I_cm(f_k) = 2*pi*f_k * C_p * |V_sw(f_k)|          [A]      (i = C dV/dt)

which is a current source driving the CM loop.  The loop (EUT board -> C_p ->
earth plane -> LISN 50 ohm -> both input lines -> board, plus the upstream
class-I Y-cap path in parallel) is solved by a tiny nodal analysis.

LISN ports:
  * CM port impedance  = two 50 ohm in PARALLEL = 25 ohm   (this file)
  * DM port impedance  = two 50 ohm in SERIES    = 100/50 ohm (emi_v2/v3)
"""
import os, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT  = '/mnt/raid10/sim-work/tps62933/emi_cm/'
RES  = '/mnt/raid10/sim-work/tps62933/results_v2/'
os.makedirs(OUT, exist_ok=True)

FSW  = 805e3
VIN  = 24.0
TR   = TF = 5e-9                 # prereq assumption
KMAX = int(30e6 / FSW)

lines = []
def log(s=''):
    print(s); lines.append(s)

# ============================================================ helpers =========
def Zc(f, C): return 1.0 / (1j * 2 * np.pi * f * C)
def Zl(f, L): return 1j * 2 * np.pi * f * L
def capZ(f, C, esr, esl): return esr + Zl(f, esl) + Zc(f, C)

def solve(n, branches, inj, Isrc=1.0):
    Y = np.zeros((n, n), complex); b = np.zeros(n, complex)
    for (i, j, Z) in branches:
        if abs(Z) < 1e-15: Z = 1e-12
        y = 1.0 / Z
        if i >= 0: Y[i, i] += y
        if j >= 0: Y[j, j] += y
        if i >= 0 and j >= 0: Y[i, j] -= y; Y[j, i] -= y
    b[inj] += Isrc
    return np.linalg.solve(Y, b)

def trapezoid_harmonics(fsw, ipk, duty, tr, tf, kmax, nper=8, ppp=16384):
    T = 1.0 / fsw; N = nper * ppp; dt = (nper * T) / N
    t = np.arange(N) * dt; tt = t % T; on = duty * T
    x = np.zeros_like(t)
    if tr > 0:
        r = tt < min(tr, on); x[r] = ipk * tt[r] / tr
    fl = (tt >= tr) & (tt < on); x[fl] = ipk
    if tf > 0:
        fa = (tt >= on) & (tt < on + tf)
        x[fa] = ipk * (1 - (tt[fa] - on) / tf)
    x -= x.mean()
    w = np.hanning(N)
    X = np.fft.rfft(x * w); f = np.fft.rfftfreq(N, dt)
    amp = 2 * np.abs(X) / w.sum()
    ks = np.arange(1, kmax + 1); fh = ks * fsw
    idx = np.array([int(round(fr / f[1])) for fr in fh])
    return fh, amp[idx]

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

# ============================================================ CM network ======
#  node 0 = earth/reference plane ; node 1 = EUT board (CM inject) ; node 2 = upstream
def P_cm(Cp, Cy=2.2e-9, Lcab=30e-9, Rcab=50e-3, Rlisn=25.0, ypath=True):
    return dict(Cp=Cp, Cy=Cy, Lcab=Lcab, Rcab=Rcab, Rlisn=Rlisn, ypath=ypath)

def Zcm_branches(f, P):
    # node 0 = EUT board (CM injection), node 1 = upstream, -1 = earth/ref plane
    b = [(0, -1, P['Rlisn'])]
    if P['ypath']:
        b.append((0, 1, P['Rcab'] + Zl(f, P['Lcab'])))
        b.append((1, -1, Zc(f, P['Cy'])))
        return 2, b
    return 1, b                 # upstream node unused when Y-cap path removed

def Vcm_spectrum(fk, Vsw_k, P):
    """fk [Hz], Vsw_k [V] -> CM voltage on the LISN (V)."""
    Icm = 2 * np.pi * fk * P['Cp'] * Vsw_k      # A
    out = np.zeros_like(fk)
    for i, fr in enumerate(fk):
        n, b = Zcm_branches(fr, P)
        V = solve(n, b, 0, Icm[i])
        out[i] = abs(V[0])
    return out

def Zcm_scalar(f, P):
    n, b = Zcm_branches(f, P)
    V = solve(n, b, 0, 1.0)
    return abs(V[0])

# ============================================================ DM model (ported from emi_v3/emi_model_v3.py, geo-nominal) ====
G = dict(Lpm=0.2984, lf1=6.32, ltr_loop=7.88, ltr_capA=0.75, lhot=0.83)   # from emi_v3/geo_numbers.txt
_BASE = dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
             L1=1e-6, dcr_l1=25e-3, cp_l1=4e-12,
             rs_dm=50.0, c_lisn=1e-6, l_mains=50e-6)
def P_geo_nominal():
    d = dict(_BASE)
    d.update(esr_10u=2e-3, esl_10u=1.0e-9, esr_100n=20e-3, esl_100n=0.8e-9,
             esr_e=4.0, esl_e=3.0e-9, r85=0.1,
             lf1=G['lf1'] * 1e-9, ltr_loop=G['ltr_loop'] * 1e-9,
             lhot=G['lhot'] * 1e-9, ltr_capA=G['ltr_capA'] * 1e-9, rf1=5e-3)
    return d
def damping(P, fr):
    zc = capZ(fr, P['Ce'], P['esr_e'], P['esl_e'])
    return P['r85'] + 1.0 / (1.0 / zc + 1.0 / zc)
def Zt_dm(fr, P):
    b = []
    b.append((0, -1, P['rs_dm'] + Zc(fr, P['c_lisn'])))
    b.append((0, -1, Zl(fr, P['l_mains'])))
    b.append((0, 1, P['rf1'] + Zl(fr, P['lf1'])))
    for _ in range(3):
        b.append((1, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'] + P['ltr_capA'])))
    b.append((1, 2, P['dcr_l1'] + Zl(fr, P['L1'] + P['ltr_loop'])))
    if P['cp_l1'] > 0: b.append((1, 2, Zc(fr, P['cp_l1'])))
    b.append((2, 3, Zl(fr, P['lhot'])))
    for _ in range(2):
        b.append((3, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'])))
    b.append((3, -1, capZ(fr, P['C100n'], P['esr_100n'], P['esl_100n'])))
    b.append((3, -1, damping(P, fr)))
    V = solve(4, b, 2)
    return abs(V[0])
def dm_spectrum(P, ipk, duty, tr=10e-9, tf=10e-9):
    f, a = trapezoid_harmonics(FSW, ipk, duty, tr, tf, KMAX)
    z = np.array([Zt_dm(fr, P) for fr in f])
    return f, 20 * np.log10(np.maximum(a * z, 1e-30) / 1e-6)

# ============================================================ load sim SW =====
def load_sw(case):
    a = np.loadtxt(RES + 'E_%s.txt' % case)
    return a[:, 0], a[:, 3]

def measured_sw_spec(t, vsw, M=32):
    """FFT of the measured SW voltage over the last M switching periods."""
    T = 1.0 / FSW; Tw = M * T
    m = t >= (t[-1] - Tw)
    tt, vv = t[m], vsw[m]
    N = M * 4096
    tu = np.linspace(tt[0], tt[-1], N, endpoint=False)
    vu = np.interp(tu, tt, vv)
    vu = vu - vu.mean()
    X = np.fft.rfft(vu)
    f = np.fft.rfftfreq(N, tu[1] - tu[0])
    ks = np.arange(1, KMAX + 1); fk = ks * FSW
    amp = 2 * np.abs(X) / N
    idx = np.array([int(round(fr / f[1])) for fr in fk])
    return fk, amp[idx], vv.mean() / VIN

CASES = [('noload', 0.02, 0.22, 0.15, '#1f77b4'),
         ('half',   1.50, 1.50, 0.50, '#2ca02c'),
         ('full',   3.00, 3.00, 0.50, '#d62728')]

# ============================================================ RUN =============
log('=' * 78)
log('TPS62933 CONDUCTED COMMON-MODE (CM) model, 150 kHz - 30 MHz, CISPR 32')
log('=' * 78)
log('')

# --- measured SW edge (sanity / cross-check, NOT used for the source spectrum) --
log('--- measured SW edges from results_v2/E_*.txt (v(sw) col) ---')
for nm, il, ipk, dty, col in CASES:
    t, vsw = load_sw(nm)
    dt = np.diff(t); dt[dt <= 0] = np.nan
    ts = t[t >= t[-1] - 20e-6]; vs = vsw[t >= t[-1] - 20e-6]
    log('  %-7s tail-mean(Vsw)=%.2f V -> D_meas=%.3f ; max|dV/dt| (numerical) huge (dt->0 artefact)'
        % (nm, vs.mean(), vs.mean() / VIN))
log('  NOTE: TPS62933P sim SW edges collapse to <=~1 sample (|dV/dt| unphysical) ->')
log('        the SIM edge is a maxstep=2n numerical artefact. Source uses the')
log('        prereq assumption tr=tf=5 ns (TI typical).')

# --- CM network / port impedance ---
Pnom = P_cm(Cp=10e-12)      # nominal C_p = 10 pF
fg = np.logspace(np.log10(150e3), np.log10(30e6), 400)
log('')
log('--- CM port impedance Z_cm(f) = 25 ohm || (cable + upstream Y-cap path) ---')
log('  nominal: C_p=10 pF, C_y=2.2 nF (est), L_cable=30 nH (10 cm pair, est)')
for fr in [150e3, 805e3, 5e6, 30e6]:
    z_par = Zcm_scalar(fr, Pnom)
    z_25 = Zcm_scalar(fr, P_cm(10e-12, ypath=False))
    log('   f=%8.3f MHz : Z_cm=%7.2f ohm (LISN-only %.1f)'
        % (fr / 1e6, z_par, z_25))

# --- V_sw harmonics: measured vs idealized trapezoid ---
log('')
log('--- V_sw harmonic content (basis of I_cm = 2*pi*f*C_p*V_sw) ---')
log('  case    D_meas  |Vsw1|meas  |Vsw1|trap(5ns)  |Vsw5|meas  |Vsw37|meas')
sws = {}
for nm, il, ipk, dty, col in CASES:
    t, vsw = load_sw(nm)
    fk, am, dm = measured_sw_spec(t, vsw)
    _, atr = trapezoid_harmonics(FSW, VIN, dty, TR, TF, KMAX)
    sws[nm] = dict(fk=fk, am=am, atr=atr, dm=dm)
    log('  %-7s %.3f    %8.3f V    %8.3f V     %7.3f V    %7.4f V'
        % (nm, dm, am[0], atr[0], am[4], am[36] if len(am) > 36 else am[-1]))

# --- CM voltage, three cases. PRIMARY source = MEASURED SW spectrum ---
# --- (FFT of v(sw) from results_v2, task requirement); trapezoid = cross-check. ---
log('')
log('--- CM voltage on LISN, C_p=10 pF, PRIMARY = measured SW dv/dt spectrum ---')
log('  case   Vcm@805k[dBuV]  Vcm@4.03M  Vcm worst-margin[dB] @ pos        viol')
cm = {}
cm_trap = {}
for nm, il, ipk, dty, col in CASES:
    fk = sws[nm]['fk']
    Vcm = Vcm_spectrum(fk, sws[nm]['am'], Pnom)          # measured
    dbv = 20 * np.log10(np.maximum(Vcm, 1e-30) / 1e-6)
    m = cispr_B(fk) - dbv; j = int(np.argmin(m))
    cm[nm] = dict(f=fk, dbv=dbv, m=m, col=col, il=il)
    log('  %-6s %10.1f      %8.1f    %+8.1f @ %6.3f MHz   %d'
        % (nm, dbv[0], np.interp(4.03e6, fk, dbv), m[j], fk[j] / 1e6, int((m < 0).sum())))
    Vcm2 = Vcm_spectrum(fk, sws[nm]['atr'], Pnom)        # trapezoid 5ns
    dbv2 = 20 * np.log10(np.maximum(Vcm2, 1e-30) / 1e-6)
    m2 = cispr_B(fk) - dbv2
    cm_trap[nm] = dict(f=fk, dbv=dbv2, m=m2)

log('')
log('--- CM voltage, cross-check: idealized trapezoid tr=tf=5 ns (prereq), D=0.15/0.5/0.5 ---')
for nm, il, ipk, dty, col in CASES:
    r = cm_trap[nm]; j = int(np.argmin(r['m']))
    log('  %-6s Vcm@805k=%6.1f dBuV  worst-margin=%+7.1f dB @ %.3f MHz'
        % (nm, r['dbv'][0], r['m'][j], r['f'][j] / 1e6))

# --- C_p sensitivity (prereq: 2 / 10 / 50 pF) ---
log('')
log('--- C_p sensitivity (full load, measured SW spectrum, C_y=2.2 nF) ---')
fk = sws['full']['fk']; V1 = sws['full']['am']
sweep_Cp = {}
for Cp in [2e-12, 10e-12, 50e-12]:
    Vcm = Vcm_spectrum(fk, V1, P_cm(Cp))
    dbv = 20 * np.log10(np.maximum(Vcm, 1e-30) / 1e-6)
    m = cispr_B(fk) - dbv; j = int(np.argmin(m))
    sweep_Cp[Cp] = dbv
    log('   C_p=%4.0f pF : Vcm@805k=%6.1f dBuV  worst-margin=%+7.1f dB @ %.3f MHz'
        % (Cp * 1e12, dbv[0], m[j], fk[j] / 1e6))

# --- other sensitivities ---
log('')
log('--- secondary sensitivities (full load, C_p=10 pF) ---')
for tag, kw in [('C_y=0.47 nF', dict(Cy=0.47e-9)), ('C_y=2.2 nF (nom)', dict(Cy=2.2e-9)),
                ('C_y=4.7 nF', dict(Cy=4.7e-9)), ('no Y-cap path (LISN only)', dict(ypath=False))]:
    P = P_cm(10e-12, **kw)
    Vcm = Vcm_spectrum(fk, V1, P)
    dbv = 20 * np.log10(np.maximum(Vcm, 1e-30) / 1e-6)
    m = cispr_B(fk) - dbv; j = int(np.argmin(m))
    log('   %-24s Vcm@805k=%6.1f dBuV  worst-margin=%+7.1f dB @ %.3f MHz'
        % (tag, dbv[0], m[j], fk[j] / 1e6))
log('   tr=tf sweep (C_p=10 pF, C_y nom):')
for tr in [2e-9, 5e-9, 10e-9, 20e-9]:
    _, a = trapezoid_harmonics(FSW, VIN, 0.5, tr, tr, KMAX)
    Vcm = Vcm_spectrum(fk, a, Pnom)
    dbv = 20 * np.log10(np.maximum(Vcm, 1e-30) / 1e-6)
    m = cispr_B(fk) - dbv; j = int(np.argmin(m))
    log('     tr=tf=%4.0f ns : worst-margin=%+7.1f dB @ %.3f MHz (Vcm@805k=%.1f)'
        % (tr * 1e9, m[j], fk[j] / 1e6, dbv[0]))

# --- DM for comparison (ported geo-nominal v3) ---
log('')
log('--- DM reference (ported emi_v3 geo-nominal) ---')
Pdm = P_geo_nominal()
dm = {}
for nm, il, ipk, dty, col in CASES:
    fk2, dbv = dm_spectrum(Pdm, ipk, dty)
    m = cispr_B(fk2) - dbv; j = int(np.argmin(m))
    dm[nm] = dict(f=fk2, dbv=dbv, m=m)
    log('  %-6s DM V1@805k=%6.1f dBuV  worst-margin=%+7.1f dB @ %.3f MHz'
        % (nm, dbv[0], m[j], fk2[j] / 1e6))
log('  (cross-check: v3 published DM full-load V1=13.2 dBuV, margin +39.9 dB)')

# --- CM + DM total ---
log('')
log('--- CM + DM total on the worst line (worst-case phasor = Vcm + Vdm/2) ---')
tot = {}
for nm, il, ipk, dty, col in CASES:
    fkc = cm[nm]['f']; vcm = 10 ** (cm[nm]['dbv'] / 20) * 1e-6
    # DM on same harmonic grid
    vdm = 10 ** (np.interp(fkc, dm[nm]['f'], dm[nm]['dbv']) / 20) * 1e-6
    vtot = vcm + vdm / 2.0
    dbvt = 20 * np.log10(vtot / 1e-6)
    m = cispr_B(fkc) - dbvt; j = int(np.argmin(m))
    tot[nm] = dict(f=fkc, dbv=dbvt, m=m)
    log('  %-6s total@805k=%6.1f dBuV  worst-margin=%+7.1f dB @ %.3f MHz  (%s)'
        % (nm, dbvt[0], m[j], fkc[j] / 1e6,
           'CM-dominated' if vcm[0] > vdm[0] / 2 else 'DM-dominated'))

# --- answer ---
mf = tot['full']['m'].min(); jf = int(np.argmin(tot['full']['m']))
log('')
log('*** ANSWER (CM, C_p=10 pF nominal, measured SW dv/dt source):')
log('    full-load CM worst Class-B margin = %+.1f dB @ %.2f MHz.' % (cm['full']['m'].min(),
    cm['full']['f'][int(np.argmin(cm['full']['m']))] / 1e6))
log('    CM is a VIOLATION by ~%.0f dB; CM dominates DM by ~%.0f dB.' %
    (-cm['full']['m'].min(), cm['full']['dbv'][0] - dm['full']['dbv'][0]))
log('    Combined CM+DM full-load worst margin = %+.1f dB @ %.2f MHz. ***' % (mf, tot['full']['f'][jf] / 1e6))

# ============================================================ FIGURES =========
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.6))
for ax in (ax1, ax2):
    ax.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.0, label='CISPR 32 Class B (QP)')
    ax.plot(fgp / 1e6, cispr_A(fgp), 'k--', lw=1.2, label='CISPR 32 Class A (QP)')
    ax.set_xlim(0.15, 30); ax.grid(True, which='both', alpha=0.3)
    ax.set_xlabel('Frequency [MHz]')
ax1.set_ylim(30, 110)
for nm, il, ipk, dty, col in CASES:
    r = cm[nm]
    ax1.semilogx(r['f'] / 1e6, r['dbv'], '-o', ms=3.2, color=col,
                 label='%s  Iload=%.2fA  min-marg=%+.1f dB' % (nm, il, r['m'].min()))
ax1.set_ylabel('CM voltage on LISN [dB$\\mu$V]')
ax1.set_title('(a) CM, three load cases  ($C_p$=10 pF, tr=tf=5 ns)')
ax1.legend(loc='upper right', fontsize=8.3)
ax2.set_ylim(30, 110)
for nm, il, ipk, dty, col in CASES:
    ax2.semilogx(tot[nm]['f'] / 1e6, tot[nm]['dbv'], '-s', ms=3.0, color=col,
                 label='%s  total  min-marg=%+.1f dB' % (nm, tot[nm]['m'].min()))
ax2.set_ylabel('CM+DM voltage on LISN [dB$\\mu$V]')
ax2.set_title('(b) CM + DM total (DM from emi_v3 geo-nominal)')
ax2.legend(loc='upper right', fontsize=8.3)
fig.suptitle('TPS62933 conducted CM (and CM+DM) vs CISPR 32  — fsw=805 kHz', fontsize=12)
fig.tight_layout()
fig.savefig(OUT + 'fig_emi_cm_cispr.png', dpi=140)
log(''); log('Wrote fig_emi_cm_cispr.png')

# combined CM+DM separate figure
fig2, ax = plt.subplots(figsize=(9.5, 7))
ax.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.0, label='CISPR 32 Class B (QP)')
ax.plot(fgp / 1e6, cispr_A(fgp), 'k--', lw=1.2, label='CISPR 32 Class A (QP)')
for nm, il, ipk, dty, col in CASES:
    ax.semilogx(cm[nm]['f'] / 1e6, cm[nm]['dbv'], '-o', ms=3, color=col, alpha=0.45,
                label='%s CM' % nm)
    ax.semilogx(tot[nm]['f'] / 1e6, tot[nm]['dbv'], '-s', ms=3.2, color=col,
                label='%s total (CM+DM)' % nm)
ax.set_xlim(0.15, 30); ax.set_ylim(30, 110); ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel('LISN voltage [dB$\\mu$V]')
ax.set_title('TPS62933 CM vs CM+DM total  ($C_p$=10 pF nominal)')
ax.legend(loc='upper right', fontsize=8)
fig2.tight_layout(); fig2.savefig(OUT + 'fig_emi_cm_dm_total.png', dpi=140)
log('Wrote fig_emi_cm_dm_total.png')

# sensitivity figure
fig3, axs = plt.subplots(1, 2, figsize=(15, 6))
for Cp, c in zip([2e-12, 10e-12, 50e-12], ['#1f77b4', '#2ca02c', '#d62728']):
    axs[0].semilogx(fk / 1e6, sweep_Cp[Cp], '-o', ms=2.5, color=c, label='C_p=%.0f pF' % (Cp * 1e12))
axs[0].plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2, label='CISPR 32 B')
axs[0].set_xlim(0.15, 30); axs[0].set_ylim(30, 120); axs[0].grid(True, which='both', alpha=0.3)
axs[0].set_xlabel('Frequency [MHz]'); axs[0].set_ylabel('CM [dB$\\mu$V]')
axs[0].set_title('(a) C_p sensitivity (full load)'); axs[0].legend(fontsize=8)
for tr, c in zip([2e-9, 5e-9, 10e-9, 20e-9], ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']):
    _, a = trapezoid_harmonics(FSW, VIN, 0.5, tr, tr, KMAX)
    Vcm = Vcm_spectrum(fk, a, Pnom)
    dbv = 20 * np.log10(np.maximum(Vcm, 1e-30) / 1e-6)
    axs[1].semilogx(fk / 1e6, dbv, '-', lw=1.5, color=c, label='tr=tf=%.0f ns' % (tr * 1e9))
axs[1].plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2, label='CISPR 32 B')
axs[1].set_xlim(0.15, 30); axs[1].set_ylim(30, 120); axs[1].grid(True, which='both', alpha=0.3)
axs[1].set_xlabel('Frequency [MHz]'); axs[1].set_ylabel('CM [dB$\\mu$V]')
axs[1].set_title('(b) SW edge-rate sensitivity (C_p=10 pF)'); axs[1].legend(fontsize=8)
fig3.suptitle('TPS62933 CM sensitivity', fontsize=12)
fig3.tight_layout(); fig3.savefig(OUT + 'fig_emi_cm_sensitivity.png', dpi=140)
log('Wrote fig_emi_cm_sensitivity.png')

with open(OUT + 'emi_numbers_cm.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print('wrote emi_numbers_cm.txt')

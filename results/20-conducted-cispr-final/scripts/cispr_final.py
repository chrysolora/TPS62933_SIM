#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cispr_final.py -- TPS62933 conducted CISPR 32 Class B, OFFICIAL combined figure
                 DM + CM on ONE figure (QP), AV on a second panel.
Writes ONLY to /mnt/raid10/sim-work/tps62933/cispr_final/.

Unified (single-provenance) bus:
  * DM  : model = dm_recheck/emi_model_v2_repro.py (READ-ONLY reuse),
          inductances set to THIS-BOARD REAL COPPER geometry:
            L_hot  = 2.70 nH (microstrip over plane, l=9.3 w=3.0 h=1.43 er4.4)
                     .. 4.44 nH (Grover strip, no-return upper form)
            L_tr_capA = 2.70 nH  (30uF-bank segment, same geometry)
          source = analytic trapezoid @ fsw=805k, tr=tf=10ns, per-case (ipk,D)
          (identical to dm_real.py / emi_numbers_v2 nominal).
  * CM  : model = cm_redo/cm_redo_model.py (READ-ONLY reuse), C_p = 2.416 pF
          (cm_audit real-copper geometry: 2-layer, core 1.43mm FR4 eps_r4.3,
           top coplanar GND, FullCopper).  source = FFT(v(sw)) from
          results_v2/E_*.txt.  V_cm = 25 ohm * 2*pi*f*C_p*|V_sw|.
Both vs CISPR 32 / EN 55032 Class B QP and AV limits.
"""
import os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')

OUT = '/mnt/raid10/sim-work/tps62933/cispr_final'
RES = '/mnt/raid10/sim-work/tps62933/results_v2'
os.makedirs(OUT, exist_ok=True)
FSW = 805e3
RCM = 25.0           # LISN CM equivalent (2x50 ohm ||)
M   = 32             # switching periods for v(sw) FFT

# -------- reuse the audited DM model (read-only) ------------------------------
SRC = '/mnt/raid10/sim-work/tps62933/dm_recheck/dm_recheck.py'
_src = open(SRC).read()
_prefix = _src.split('lines = []')[0]
G = {}
exec(_prefix, G)
mk = G['mk']; P_nominal = G['P_nominal']
worst_margin = G['worst_margin']; cispr_B = G['cispr_B']; cispr_B_avg = G['cispr_B_avg']
spectrum = G['spectrum']

# real copper geometry inductance (dm_real.py formulas, verbatim)
def grover(l, w, t=0.03):
    return 0.2 * l * (np.log(2*l/(w+t)) + 0.5 + 0.2235*(w+t)/l)
def microstrip(l, w, h=1.43, er=4.4):
    w_h = w/h
    eeff = (er+1)/2 + (er-1)/2/np.sqrt(1+12/w_h)
    if w_h <= 1:
        Z0 = 60/np.sqrt(eeff)*np.log(8*h/w + w/(4*h))
    else:
        Z0 = 120*np.pi/(np.sqrt(eeff)*(w/h + 1.393 + 0.667*np.log(w/h + 1.444)))
    return Z0*np.sqrt(eeff)/2.99792458e8*1e6*l

L_HOT_MIC = microstrip(9.3, 3.0)      # 2.70 nH  (real, return-plane)
L_HOT_GRO = grover(9.3, 3.0)          # 4.44 nH  (real, no-return upper)
L_TR_MIC  = microstrip(4.0, 3.0)      # 30uF-bank seg, same geometry (<= hot)
L_TR_REAL = 2.70                      # dm_real nominal for bank seg

# -------- CM model functions (verbatim from cm_redo_model.py) -----------------
def load_sw(case):
    a = np.loadtxt(RES + '/E_%s.txt' % case)
    return a[:, 0], a[:, 3]          # t, v(sw)
def sw_harmonics(t, vsw, M=M):
    T = 1.0/FSW; Tw = M*T
    m = t >= (t[-1] - Tw)
    tt, vv = t[m], vsw[m]
    N = M*4096
    tu = np.linspace(tt[0], tt[-1], N, endpoint=False)
    vu = np.interp(tu, tt, vv); vu = vu - vu.mean()
    X = np.fft.rfft(vu); f = np.fft.rfftfreq(N, tu[1]-tu[0])
    KMAX = int(30e6/FSW)
    ks = np.arange(1, KMAX+1); fk = ks*FSW
    amp = 2*np.abs(X)/N
    idx = np.array([int(round(fr/f[1])) for fr in fk])
    return fk, amp[idx]
def dbuv(v): return 20*np.log10(np.maximum(v, 1e-30)/1e-6)

CP_NOM = 2.416e-12   # fullcopper eps_r 4.3
CP46   = 2.572e-12

CASES = [('noload', 0.02, 0.22, 0.15, '#1f77b4'),
         ('half',   1.50, 1.50, 0.50, '#2ca02c'),
         ('full',   3.00, 3.00, 0.50, '#d62728')]

L = []
def log(s=''):
    print(s); L.append(s)

log('='*78)
log('TPS62933 CONDUCTED CISPR 32 Class B -- UNIFIED DM + CM (official)')
log('='*78)
log('DM model : dm_recheck/emi_model_v2_repro.py (read-only reuse)')
log('CM model : cm_redo/cm_redo_model.py  (read-only reuse)')
log('Real-copper L_hot (24VIN 3.0x9.3mm over GND plane 1.43mm):')
log('   microstrip-over-plane = %.2f nH ; Grover strip = %.2f nH  -> 2.7-4.4 nH'
    % (L_HOT_MIC, L_HOT_GRO))
log('')

# ---------------- DM per case: real copper geometry ---------------------------
log('--- DM (REAL copper geometry) : 3 load cases, Class B ---')
log('%-7s %5s %5s  %8s %8s  %-22s  %-22s'
    % ('case', 'Ild', 'D', 'Ipk[A]', '|Zt|1', 'QP margin', 'AV margin'))
dm = {}
for nm, il, ipk, duty, col in CASES:
    P = mk(P_nominal()); P['lhot'] = L_HOT_MIC*1e-9; P['ltr_capA'] = L_TR_REAL*1e-9
    f, a, z, dbv = spectrum('pi', P, ipk, duty, 10e-9, 10e-9)
    mq = cispr_B(f) - dbv; ma = cispr_B_avg(f) - dbv
    jq = int(np.argmin(mq)); ja = int(np.argmin(ma))
    dm[nm] = dict(f=f, dbv=dbv, mq=mq, ma=ma, col=col, il=il)
    log('%-7s %5.2f %5.2f  %8.3f %8.1e  QP %+6.1f @%5.2fMHz  AV %+6.1f @%5.2fMHz'
        % (nm, il, duty, ipk, z[0], mq[jq], f[jq]/1e6, ma[ja], f[ja]/1e6))
# DM with Grover (conservative real) hot-loop L -> show range
for nm, il, ipk, duty, col in CASES:
    P = mk(P_nominal()); P['lhot'] = L_HOT_GRO*1e-9; P['ltr_capA'] = L_TR_REAL*1e-9
    m, fq, _, _ = worst_margin('pi', P, ipk, duty, 10e-9, 10e-9, lim=cispr_B)
    ma_, _, _, _ = worst_margin('pi', P, ipk, duty, 10e-9, 10e-9, lim=cispr_B_avg)
    dm[nm]['qp_gro'] = m; dm[nm]['av_gro'] = ma_
    log('   [L_hot=%.2fnH Grover] %-6s QP %+6.1f @%5.2fMHz  AV %+6.1f dB'
        % (L_HOT_GRO, nm, m, fq/1e6, ma_))
log('')

# ---------------- CM per case: real-copper C_p --------------------------------
log('--- CM (C_p=%.3f pF real copper, primary topology) : 3 load cases ---' % (CP_NOM*1e12))
log('%-7s  %8s  %8s  %-22s  %-22s' % ('case', 'Vsw1[V]', 'Vcm1', 'QP margin', 'AV margin'))
cm = {}
for nm, il, ipk, duty, col in CASES:
    t, v = load_sw(nm)
    fk, am = sw_harmonics(t, v)
    d = dbuv(RCM*(2*np.pi*fk*CP_NOM*am))
    mq = cispr_B(fk) - d; ma = cispr_B_avg(fk) - d
    jq = int(np.argmin(mq)); ja = int(np.argmin(ma))
    cm[nm] = dict(f=fk, dbv=d, mq=mq, ma=ma, col=col, il=il, steady=(nm != 'half'))
    log('%-7s  %8.3f  %8.1f  QP %+6.1f @%5.2fMHz  AV %+6.1f @%5.2fMHz %s'
        % (nm, am[0], d[0], mq[jq], fk[jq]/1e6, ma[ja], fk[ja]/1e6,
           '' if nm != 'half' else '(half NOT steady-state)'))
log('')

# ---------------- numbers file ------------------------------------------------
def fmt_margin(x): return '%+.1f' % x
with open(OUT + '/numbers_cispr_final.txt', 'w') as fh:
    fh.write('\n'.join(L) + '\n')

# ---------------- FIGURE: QP (DM+CM same panel) + AV --------------------------
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)
fig, (axQ, axA) = plt.subplots(1, 2, figsize=(17, 7.2))
for ax, lim, ttl in ((axQ, cispr_B, 'QUASI-PEAK'), (axA, cispr_B_avg, 'AVERAGE')):
    ax.plot(fgp/1e6, lim(fgp), 'k-', lw=2.6, label='CISPR 32 Class B %s limit' % ttl)
    ax.set_xscale('log'); ax.set_xlim(0.15, 30); ax.set_ylim(20, 105)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel('LISN voltage [dB$\\mu$V]')
    ax.set_title('%s detector' % ttl)
# DM : solid, CM : dashed ; colours by load case
for nm, il, ipk, duty, col in CASES:
    r = dm[nm]
    axQ.semilogx(r['f']/1e6, r['dbv'], '-', color=col, lw=1.7,
                 label='DM %-6s Ild=%.1fA  QPm %+.1f dB' % (nm, il, r['mq'].min()))
    axA.semilogx(r['f']/1e6, r['dbv'], '-', color=col, lw=1.7,
                 label='DM %-6s Ild=%.1fA  AVm %+.1f dB' % (nm, il, r['ma'].min()))
for nm, il, ipk, duty, col in CASES:
    r = cm[nm]
    axQ.semilogx(r['f']/1e6, r['dbv'], '--', color=col, lw=1.7,
                 label='CM %-6s Ild=%.1fA  QPm %+.1f dB' % (nm, il, r['mq'].min()))
    axA.semilogx(r['f']/1e6, r['dbv'], '--', color=col, lw=1.7,
                 label='CM %-6s Ild=%.1fA  AVm %+.1f dB' % (nm, il, r['ma'].min()))
axQ.text(0.16, 100, 'solid = DM   dashed = CM', fontsize=8.5)
axQ.legend(loc='lower left', fontsize=7.2, ncol=2)
axA.legend(loc='lower left', fontsize=7.2, ncol=2)
fig.suptitle('TPS62933 conducted emission vs CISPR 32 Class B -- DM & CM, real-copper bus\n'
             'fsw=805 kHz | DM: L_hot=2.70 nH microstrip real | CM: C_p=2.416 pF real copper',
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT + '/fig_conducted_cispr_final.png', dpi=150)
log('Wrote fig_conducted_cispr_final.png')

with open(OUT + '/numbers_cispr_final.txt', 'w') as fh:
    fh.write('\n'.join(L) + '\n')
print('DONE')

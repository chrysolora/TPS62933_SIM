#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_cispr_split.py -- TPS62933 conducted CISPR 32 Class B, split into 3 single-panel figures.
NO re-simulation: reuses the audited DM model (dm_recheck.py prefix, read-only) and the
CM model (verbatim from cm_redo/cm_redo_model.py). Figures are drawn from those models only.
Writes ONLY to /mnt/raid10/sim-work/tps62933/cispr_split/.
"""
import os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')

OUT = '/mnt/raid10/sim-work/tps62933/cispr_split'
RES = '/mnt/raid10/sim-work/tps62933/results_v2'
os.makedirs(OUT, exist_ok=True)
FSW = 805e3
RCM = 25.0
M   = 32

# -------- reuse the audited DM model (read-only) ------------------------------
SRC = '/mnt/raid10/sim-work/tps62933/dm_recheck/dm_recheck.py'
_src = open(SRC).read()
_prefix = _src.split('lines = []')[0]
G = {}
exec(_prefix, G)
mk = G['mk']; P_nominal = G['P_nominal']
cispr_B = G['cispr_B']; cispr_B_avg = G['cispr_B_avg']
trapezoid_harmonics = G['trapezoid_harmonics']
Zc = G['Zc']; Zl = G['Zl']; capZ = G['capZ']; solve = G['solve']
damping = G['damping']; KMAX = G['KMAX']

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
L_HOT_MIC = microstrip(9.3, 3.0)      # 2.70 nH real copper (return-plane)
L_HOT_GRO = grover(9.3, 3.0)          # 4.44 nH real copper (no-return upper)
L_TR_REAL = 2.70

def Zt_mode(net, fr, P):
    b = []
    b.append((0, -1, P['rs_dm'] + Zc(fr, P['c_lisn'])))
    b.append((0, -1, Zl(fr, P['l_mains'])))
    b.append((0, 1, P['rf1'] + Zl(fr, P['lf1'])))
    for _ in range(3):
        b.append((1, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'] + P['ltr_capA'])))
    if net == 'pi':
        L1, dcr, cpl1, lloop = P['L1'], P['dcr_l1'], P['cp_l1'], P['ltr_loop']
    else:
        L1, dcr, cpl1, lloop = 0.0, 0.0, 0.0, 0.0
    b.append((1, 2, dcr + Zl(fr, L1 + lloop)))
    if cpl1 > 0: b.append((1, 2, Zc(fr, cpl1)))
    b.append((2, 3, Zl(fr, P['lhot'])))
    if net != 'nopi_b':
        for _ in range(2):
            b.append((3, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'])))
        b.append((3, -1, capZ(fr, P['C100n'], P['esr_100n'], P['esl_100n'])))
        b.append((3, -1, damping(P, fr)))
    V = solve(4, b, 2)
    return abs(V[0])

def spectrum_mode(net, P, ipk, duty, tr, tf, kmax=KMAX):
    f, a = trapezoid_harmonics(FSW, ipk, duty, tr, tf, kmax)
    z = np.array([Zt_mode(net, fr, P) for fr in f])
    return f, a, z, 20*np.log10(np.maximum(a*z, 1e-30)/1e-6)

def dm_case(net, ipk, duty=0.5):
    P = mk(P_nominal()); P['lhot'] = L_HOT_MIC*1e-9; P['ltr_capA'] = L_TR_REAL*1e-9
    f, a, z, dbv = spectrum_mode(net, P, ipk, duty, 10e-9, 10e-9)
    mq = cispr_B(f) - dbv; ma = cispr_B_avg(f) - dbv
    return dict(f=f, dbv=dbv, mq=mq, ma=ma, jq=int(np.argmin(mq)), ja=int(np.argmin(ma)))

# -------- CM model (verbatim from cm_redo_model.py) ---------------------------
def load_sw(case):
    a = np.loadtxt(RES + '/E_%s.txt' % case)
    return a[:, 0], a[:, 3]
def sw_harmonics(t, vsw, M=M):
    T = 1.0/FSW; Tw = M*T
    m = t >= (t[-1] - Tw)
    tt, vv = t[m], vsw[m]
    N = M*4096
    tu = np.linspace(tt[0], tt[-1], N, endpoint=False)
    vu = np.interp(tu, tt, vv); vu = vu - vu.mean()
    X = np.fft.rfft(vu); f = np.fft.rfftfreq(N, tu[1]-tu[0])
    KM = int(30e6/FSW)
    ks = np.arange(1, KM+1); fk = ks*FSW
    amp = 2*np.abs(X)/N
    idx = np.array([int(round(fr/f[1])) for fr in fk])
    return fk, amp[idx]
def dbuv(v): return 20*np.log10(np.maximum(v, 1e-30)/1e-6)
CP_NOM = 2.416e-12

def cm_case(case):
    t, v = load_sw(case)
    fk, am = sw_harmonics(t, v)
    d = dbuv(RCM*(2*np.pi*fk*CP_NOM*am))
    mq = cispr_B(fk) - d; ma = cispr_B_avg(fk) - d
    return dict(f=fk, dbv=d, mq=mq, ma=ma, jq=int(np.argmin(mq)), ja=int(np.argmin(ma)))

# -------- compute (full-load case = representative worst) ---------------------
DM  = {k: dm_case(k, 3.0) for k in ('pi', 'nopi_a', 'nopi_b')}
CM  = {c: cm_case(c) for c in ('noload', 'half', 'full')}

lg = []
def log(s=''):
    print(s); lg.append(s)

log('='*74)
log('TPS62933 CISPR 32 Class B -- SPLIT single-panel figures (no re-simulation)')
log('DM source : dm_recheck/dm_recheck.py prefix (read-only) | L_hot=%.2f nH real copper | full load Ipk=3A' % L_HOT_MIC)
log('CM source : cm_redo/cm_redo_model.py (verbatim) | C_p=%.3f pF real copper | v(sw) from results_v2/E_*.txt' % (CP_NOM*1e12))
log('='*74)
for k, lab in (('pi','with pi'), ('nopi_a','no pi (a) L1 shorted'), ('nopi_b','no pi (b) whole pi removed')):
    r = DM[k]
    log('DM full %-26s QP %+6.1f @%5.2fMHz   AV %+6.1f @%5.2fMHz'
        % (lab, r['mq'][r['jq']], r['f'][r['jq']]/1e6, r['ma'][r['ja']], r['f'][r['ja']]/1e6))
for c in ('noload','half','full'):
    r = CM[c]
    log('CM %-8s QP %+6.1f @%5.2fMHz   AV %+6.1f @%5.2fMHz'
        % (c, r['mq'][r['jq']], r['f'][r['jq']]/1e6, r['ma'][r['ja']], r['f'][r['ja']]/1e6))
sup_a = DM['nopi_a']['dbv'][DM['nopi_a']['jq']] - DM['pi']['dbv'][DM['pi']['jq']]
sup_b = DM['nopi_b']['dbv'][DM['nopi_b']['jq']] - DM['pi']['dbv'][DM['pi']['jq']]
log('pi-filter DM suppression (full load, at 4.03 MHz worst point):')
log('   emission(nopi_a) - emission(pi) = %+.1f dB  -> margin improves by %+.1f dB' % (sup_a, sup_a))
log('   emission(nopi_b) - emission(pi) = %+.1f dB  -> margin improves by %+.1f dB' % (sup_b, sup_b))

# -------- plotting helpers ----------------------------------------------------
FSZ = 15
plt.rcParams.update({'font.size': FSZ, 'axes.titlesize': FSZ+2, 'axes.labelsize': FSZ,
                     'legend.fontsize': FSZ-2, 'xtick.labelsize': FSZ-1, 'ytick.labelsize': FSZ-1})

def setup(ax, title):
    ax.set_xscale('log'); ax.set_xlim(0.15, 30); ax.set_ylim(0, 110)
    ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel(r'LISN voltage [dB$\mu$V]')
    ax.set_title(title)
    ax.grid(True, which='both', alpha=0.35)
    ax.set_xticks([0.15,0.5,1,2,5,10,20,30])
    ax.set_xticklabels(['0.15','0.5','1','2','5','10','20','30'])

def add_limits(ax, fg):
    ax.semilogx(fg/1e6, cispr_B(fg), 'k-', lw=3.0, label='CISPR 32 Class B QP limit')
    ax.semilogx(fg/1e6, cispr_B_avg(fg), 'k--', lw=2.2, label='CISPR 32 Class B AV limit')

def mark_worst(ax, r, color):
    fq = r['f'][r['jq']]; vq = r['dbv'][r['jq']]; mq = r['mq'][r['jq']]
    ax.plot([fq/1e6], [vq], marker='o', ms=11, mfc='none', mec='red', mew=2.5, zorder=6)
    ax.annotate('worst QP margin %+.1f dB\n@ %.2f MHz' % (mq, fq/1e6),
                xy=(fq/1e6, vq), xytext=(fq/1e6*1.35, min(vq+22, 100)),
                fontsize=FSZ-1, color='red',
                arrowprops=dict(arrowstyle='->', color='red', lw=1.6), zorder=7)

fg = np.logspace(np.log10(1.5e5), np.log10(3e7), 900)

# 1) DM no pi
fig, ax = plt.subplots(figsize=(10.5, 6.8))
r = DM['nopi_b']
ax.semilogx(r['f']/1e6, r['dbv'], color='#d62728', lw=2.6, marker='o', ms=3.2,
            label='DM noise, NO pi filter')
add_limits(ax, fg); setup(ax, 'TPS62933 DM conducted emission -- NO pi filter (full load 3.0 A)')
mark_worst(ax, r, '#d62728')
ax.legend(loc='upper right', framealpha=0.95)
fig.tight_layout(); fig.savefig(OUT + '/cispr_dm_no_pi.png', dpi=150); plt.close(fig)

# 2) DM with pi
fig, ax = plt.subplots(figsize=(10.5, 6.8))
r = DM['pi']
ax.semilogx(r['f']/1e6, r['dbv'], color='#1f77b4', lw=2.6, marker='o', ms=3.2,
            label='DM noise, WITH pi filter (L1 1uH + 20uF+100nF)')
add_limits(ax, fg); setup(ax, 'TPS62933 DM conducted emission -- WITH pi filter (full load 3.0 A)')
mark_worst(ax, r, '#1f77b4')
ax.legend(loc='upper right', framealpha=0.95)
fig.tight_layout(); fig.savefig(OUT + '/cispr_dm_with_pi.png', dpi=150); plt.close(fig)

# 3) CM (full load)
fig, ax = plt.subplots(figsize=(10.5, 6.8))
r = CM['full']
ax.semilogx(r['f']/1e6, r['dbv'], color='#2ca02c', lw=2.6, marker='o', ms=3.2,
            label='CM noise (full load 3.0 A), C_p=2.416 pF')
add_limits(ax, fg); setup(ax, 'TPS62933 CM conducted emission -- full load 3.0 A')
mark_worst(ax, r, '#2ca02c')
ax.legend(loc='upper right', framealpha=0.95)
fig.tight_layout(); fig.savefig(OUT + '/cispr_cm.png', dpi=150); plt.close(fig)

with open(OUT + '/numbers_cispr_split.txt', 'w') as fh:
    fh.write('\n'.join(lg) + '\n')
print('DONE')

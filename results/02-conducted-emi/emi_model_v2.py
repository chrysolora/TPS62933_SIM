#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
emi_model_v2.py -- TPS62933 buck : DM conducted-EMI, VERSION 2
                  = parasitics (ESR/ESL/trace-L) + realistic source + sensitivity + pi-filter value

Run:  cd /mnt/raid10/sim-work/tps62933/emi_v2 && python3 emi_model_v2.py
Out :  fig_emi_cispr_v2.png   (3 cases + ideal-vs-parasitic)
       fig_sensitivity.png    (worst-margin vs parasitic parameter)
       fig_emi_cispr_nofilter.png (pi-filter present vs removed)
       emi_numbers_v2.txt
       REPORT_emi_v2.md

Method: linear nodal analysis in the frequency domain.
  DM source = pulsating input (high-side switch) current at U21.VIN, trapezoidal
  train at fsw.  Filter + LISN solved as a small linear network with series RLC
  capacitors and explicit PCB trace inductances.  Zt(f) = |V_LISN / I_src|.

Every parasitic carries a provenance tag (DS=datasheet / TYP=typical / EST=estimate).
See REPORT_emi_v2.md for the table and the two limits of the envelope.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/mnt/raid10/sim-work/tps62933/emi_v2/'
RES = '/mnt/raid10/sim-work/tps62933/results_v2/'
FSW = 805e3
VIN, VOUT = 24.0, 12.0

def Zc(f, C): return 1.0 / (1j * 2 * np.pi * f * C)
def Zl(f, L): return 1j * 2 * np.pi * f * L
def capZ(f, C, esr, esl): return esr + Zl(f, esl) + Zc(f, C)

# ---------------------------------------------------------------- source -----
def trapezoid_harmonics(fsw, ipk, duty, tr, tf, kmax, nper=8, ppp=16384):
    """harmonic amplitudes (A) of a DC-removed trapezoidal pulse train."""
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

# ------------------------------------------------------------ network -------
def solve(n, branches, inj):
    """nodal solve; nodes indexed 0..n-1 (=GND..), inj 0-based current node."""
    Y = np.zeros((n, n), complex); b = np.zeros(n, complex)
    for (i, j, Z) in branches:
        if abs(Z) < 1e-15: Z = 1e-12        # ideal short floor
        y = 1.0 / Z
        if i >= 0: Y[i, i] += y
        if j >= 0: Y[j, j] += y
        if i >= 0 and j >= 0: Y[i, j] -= y; Y[j, i] -= y
    b[inj] += 1.0
    return np.linalg.solve(Y, b)

def damping(P, fr):
    zc = capZ(fr, P['Ce'], P['esr_e'], P['esl_e'])
    return P['r85'] + 1.0 / (1.0 / zc + 1.0 / zc)     # R85 serial C16||C64

def Zt(net, fr, P):
    """returns |V_LISN| for 1 A injected at the U21.VIN source node."""
    b = []
    b.append((0, -1, P['rs_dm'] + Zc(fr, P['c_lisn'])))   # LISN 50dm + 1uF  -> GND
    b.append((0, -1, Zl(fr, P['l_mains'])))               # LISN 50uH      -> mains/GND
    if net == 'pi':
        # nodes: 0=A(LISN/PPHV), 1=P(30uF bank), 2=V(VIN pin=source), 3=B(local caps)
        b.append((0, 1, P['rf1'] + Zl(fr, P['lf1'])))                 # F1/trace A->P
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
    else:                                                 # 'nopi'
        # nodes: 0=A(LISN/24VIN), 1=V(VIN pin=source), 2=B(local caps)
        b.append((0, 1, P['rf1'] + Zl(fr, P['lf1'])))
        b.append((1, 2, Zl(fr, P['lhot'])))
        for _ in range(2):
            b.append((2, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'])))
        b.append((2, -1, capZ(fr, P['C100n'], P['esr_100n'], P['esl_100n'])))
        b.append((2, -1, damping(P, fr)))
        V = solve(3, b, 1)
    return abs(V[0])

def Zt_vec(net, f, P):
    return np.array([Zt(net, fr, P) for fr in np.atleast_1d(f)])

# ------------------------------------------------------------ limits ---------
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

# ------------------------------------------------------------ parameter sets -
def P_v1repro():
    """exactly the v1 model: caps with ESR/ESL, NO trace L, single 50ohm LISN, tr=15n."""
    return dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
                esr_10u=3e-3, esl_10u=1.2e-9, esr_100n=20e-3, esl_100n=1.0e-9,
                esr_e=0.0, esl_e=0.0, r85=0.1,
                L1=1e-6, dcr_l1=30e-3, cp_l1=1e-12,
                lf1=0.0, ltr_loop=0.0, lhot=0.0, ltr_capA=0.0, rf1=0.0,
                rs_dm=50.0, c_lisn=1e-6, l_mains=50e-6)

def P_ideal():
    """textbook ideal: no ESR/ESL anywhere, no trace L, 50ohm LISN (optimistic bound)."""
    return dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
                esr_10u=0.0, esl_10u=0.0, esr_100n=0.0, esl_100n=0.0,
                esr_e=0.0, esl_e=0.0, r85=0.1,
                L1=1e-6, dcr_l1=0.0, cp_l1=0.0,
                lf1=0.0, ltr_loop=0.0, lhot=0.0, ltr_capA=0.0, rf1=0.0,
                rs_dm=50.0, c_lisn=1e-6, l_mains=50e-6)

def P_nominal():
    """central estimates WITH estimated PCB trace inductance."""
    return dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
                esr_10u=2e-3, esl_10u=1.0e-9, esr_100n=20e-3, esl_100n=0.8e-9,
                esr_e=4.0, esl_e=3.0e-9, r85=0.1,
                L1=1e-6, dcr_l1=25e-3, cp_l1=4e-12,
                lf1=5e-9, ltr_loop=15e-9, lhot=8e-9, ltr_capA=3e-9, rf1=5e-3,
                rs_dm=50.0, c_lisn=1e-6, l_mains=50e-6)

def P_pessimistic():
    """worst credible case: high ESL + long traces + strict 100ohm DM LISN + slow edge."""
    return dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
                esr_10u=6e-3, esl_10u=1.5e-9, esr_100n=40e-3, esl_100n=1.2e-9,
                esr_e=8.0, esl_e=5.0e-9, r85=0.1,
                L1=1e-6, dcr_l1=30e-3, cp_l1=6e-12,
                lf1=10e-9, ltr_loop=30e-9, lhot=20e-9, ltr_capA=8e-9, rf1=10e-3,
                rs_dm=100.0, c_lisn=1e-6, l_mains=50e-6)

def mk(p):
    P = dict(p)
    P['C10u'] = 10e-6; P['C100n'] = 100e-9; P['Ce'] = 10e-6
    return P

# ------------------------------------------------------------ cases ----------
# (name, Iload[A], Ipk[A] ~= Iload in CCM, duty, colour)
CASES = [('noload', 0.02, 0.22, 0.15, '#1f77b4'),
         ('half',   1.50, 1.50, 0.50, '#2ca02c'),
         ('full',   3.00, 3.00, 0.50, '#d62728')]
KMAX = int(30e6 / FSW)

def spectrum(net, P, ipk, duty, tr, tf):
    f, a = trapezoid_harmonics(FSW, ipk, duty, tr, tf, KMAX)
    z = Zt_vec(net, f, P)
    return f, a, z, 20 * np.log10(np.maximum(a * z, 1e-30) / 1e-6), 20 * np.log10(np.maximum(z, 1e-30) / 1e-6)

def worst_margin(net, P, ipk, duty, tr, tf, lim=cispr_B):
    f, a, z, dbv, _ = spectrum(net, P, ipk, duty, tr, tf)
    m = lim(f) - dbv
    j = int(np.argmin(m))
    return m[j], f[j], dbv[j], lim(f)[j], dbv, f, m

lines = []
def log(s=''):
    print(s); lines.append(s)

# ============================================================ RUN ============
log('=' * 78)
log('TPS62933 DM conducted-EMI  model v2  (fsw = %.1f kHz, parasitics + real source)' % (FSW / 1e3))
log('=' * 78)

TR_NOM, TF_NOM = 10e-9, 10e-9

# --- validation: v1-repro transimpedance ---
Pv1 = mk(P_v1repro())
zt1 = Zt('pi', 805e3, Pv1)
log('')
log('[validation] v1-repro |Zt|@805kHz = %.3g ohm  (v1 published: 6.33e-06)' % zt1)
log('             v1-repro |Zt|@30MHz  = %.3g ohm  (v1 published: 1.97e-05)' % Zt('pi', 30e6, Pv1))

# --- nominal run of 3 cases ---
Pn = mk(P_nominal())
log('')
log('--- nominal (central estimates, pi filter, 50ohm LISN, tr=tf=10ns) ---')
log('%-7s %-6s %-5s %8s %11s %10s %9s' % ('case', 'Iload', 'D', 'I_src1[A]', '|Zt|1[ohm]', 'V1[dBuV]', 'marg1[dB]'))
nom = {}
for name, il, ipk, duty, col in CASES:
    f, a, z, dbv, zdbv = spectrum('pi', Pn, ipk, duty, TR_NOM, TF_NOM)
    m = cispr_B(f) - dbv
    nom[name] = dict(f=f, a=a, z=z, dbv=dbv, zdbv=zdbv, m=m, col=col, il=il, duty=duty)
    log('%-7s %-6.2f %-5.2f %8.3f %11.3g %10.1f %9.1f' % (name, il, duty, a[0], z[0], dbv[0], m[0]))

log('')
log('Worst-case margin vs CISPR 32 Class B (min over 150k-30M), nominal:')
for name, il, ipk, duty, col in CASES:
    r = nom[name]; j = int(np.argmin(r['m']))
    log('  %-7s margin_min = %6.1f dB @ %.2f MHz (V=%.1f dBuV, lim=%.1f)  viol=%d'
        % (name, r['m'][j], r['f'][j] / 1e6, r['dbv'][j], cispr_B(r['f'])[j], int((r['m'] < 0).sum())))

# --- envelope: pessimistic (full load) ---
Pp = mk(P_pessimistic())
log('')
log('--- ENVELOPE, full load (Ipk=3A, D=0.5) ---')
env = {}
for tag, P, tr in [('ideal', mk(P_ideal()), 15e-9), ('v1-repro', Pv1, 15e-9),
                   ('nominal', Pn, 10e-9), ('pessimistic', Pp, 20e-9)]:
    m, fj, vj, limj, dbv, f, mall = worst_margin('pi', P, 3.0, 0.50, tr, tr)
    env[tag] = dict(m=m, f=fj, v=vj, dbv=dbv, f_all=f, lim=limj)
    log('  %-12s worst margin = %6.1f dB @ %6.2f MHz (V=%6.1f dBuV, lim=%4.1f)  |Zt|@805k=%.3g ohm'
        % (tag, m, fj / 1e6, vj, limj, Zt('pi', 805e3, P)))
opt_m = max(env[t]['m'] for t in env); pes_m = min(env[t]['m'] for t in env)
log('  => ENVELOPE (best=ideal/v1 .. worst=pessimistic): %+.1f dB .. %+.1f dB  (spread %.1f dB)'
    % (opt_m, pes_m, opt_m - pes_m))

# --- dominant-parasitic attribution (full load) ---
log('')
log('--- single-parasitic attribution (full load; each item -> pessimistic value alone) ---')
def clone(d): return dict(d)
attrib = [('esl 10uF bank',    'esl_10u', 1.5e-9),
          ('trace loop (L1 br)','ltr_loop', 30e-9),
          ('hot-loop trace',    'lhot', 20e-9),
          ('F1/input trace',    'lf1', 10e-9),
          ('30uF bank trace',   'ltr_capA', 8e-9),
          ('elec ESR',          'esr_e', 8.0),
          ('LISN 50->100 ohm',  'rs_dm', 100.0)]
for lab, key, val in attrib:
    Pt = clone(Pn); Pt[key] = val
    m, fj, vj, limj, *_ = worst_margin('pi', Pt, 3.0, 0.50, TR_NOM, TF_NOM)
    log('  %-20s worst margin %6.1f dB @ %6.2f MHz   (delta vs nom %+.1f dB)'
        % (lab, m, fj / 1e6, m - env['nominal']['m']))

# --- pi-filter value: no-pi comparison ---
Pnopi = mk(P_nominal())
log('')
log('--- pi-filter value (full load): with vs without pi filter (3x10uF bank + L1 removed) ---')
pi_m, pi_f, pi_v, *_ = worst_margin('pi', Pn, 3.0, 0.50, TR_NOM, TF_NOM)
nf_m, nf_f, nf_v, *_ = worst_margin('nopi', Pnopi, 3.0, 0.50, TR_NOM, TF_NOM)
log('  WITH pi    : |Zt|@805k=%.3g ohm, V1=%.1f dBuV, worst margin %+.1f dB @ %.2f MHz'
    % (Zt('pi', 805e3, Pn), 20 * np.log10(max(trapezoid_harmonics(FSW, 3.0, .5, TR_NOM, TF_NOM, 1)[1][0] * Zt('pi', 805e3, Pn), 1e-30) / 1e-6), pi_m, pi_f / 1e6))
log('  WITHOUT pi : |Zt|@805k=%.3g ohm, V1=%.1f dBuV, worst margin %+.1f dB @ %.2f MHz'
    % (Zt('nopi', 805e3, Pnopi), 20 * np.log10(max(trapezoid_harmonics(FSW, 3.0, .5, TR_NOM, TF_NOM, 1)[1][0] * Zt('nopi', 805e3, Pnopi), 1e-30) / 1e-6), nf_m, nf_f / 1e6))
log('  => pi filter gain = %.1f dB (Zt ratio %.1f dB)' % (nf_m - pi_m, 20 * np.log10(Zt('nopi', 805e3, Pnopi) / Zt('pi', 805e3, Pn))))

# --- |Zt| checkpoints, nominal ---
fg = np.logspace(np.log10(150e3), np.log10(30e6), 3000)
ztg = Zt_vec('pi', fg, Pn)
log('')
log('|Zt| nominal: @150kHz=%.3g  @805kHz=%.3g  @5MHz=%.3g  @30MHz=%.3g ohm'
    % (np.interp(150e3, fg, ztg), np.interp(805e3, fg, ztg),
       np.interp(5e6, fg, ztg), np.interp(30e6, fg, ztg)))

# ============================================================ FIGS ===========
# FIG 1 : 3 cases + ideal-vs-parasitic
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.6))
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)
for ax in (ax1, ax2):
    ax.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.0, label='CISPR 32 Class B')
    ax.plot(fgp / 1e6, cispr_A(fgp), 'k--', lw=1.2, label='CISPR 32 Class A')
    ax.set_xlim(0.15, 30); ax.set_ylim(-30, 100)
    ax.grid(True, which='both', alpha=0.3); ax.set_xlabel('Frequency [MHz]')
for name, il, ipk, duty, col in CASES:
    r = nom[name]
    ax1.semilogx(r['f'] / 1e6, r['dbv'], '-o', ms=3.5, color=col,
                 label='%s Iload=%.2fA  min-marg=%+.1f dB' % (name, il, r['m'].min()))
ax1.set_ylabel('DM noise on LISN [dB$\\mu$V]')
ax1.set_title('(a) three load cases, nominal parasitics')
ax1.legend(loc='upper right', fontsize=8.5)
# ideal vs parasitic (full)
_, ai, _, dbvi, _ = spectrum('pi', mk(P_ideal()), 3.0, 0.5, 15e-9, 15e-9)
_, an, _, dbvn, _ = spectrum('pi', Pn, 3.0, 0.5, 10e-9, 10e-9)
_, ap, _, dbvp, _ = spectrum('pi', Pp, 3.0, 0.5, 20e-9, 20e-9)
ax2.semilogx(nom['full']['f'] / 1e6, dbvi, '-o', ms=3, color='#2ca02c', label='ideal (no parasitics, no trace L)')
ax2.semilogx(nom['full']['f'] / 1e6, dbvn, '-s', ms=3, color='#1f77b4', label='nominal parasitics')
ax2.semilogx(nom['full']['f'] / 1e6, dbvp, '-^', ms=3, color='#d62728', label='pessimistic parasitics')
ax2.set_ylabel('DM noise on LISN [dB$\\mu$V]')
ax2.set_title('(b) full load: parasitics envelope')
ax2.legend(loc='upper right', fontsize=8.5)
fig.suptitle('TPS62933 DM conducted emission v2 (fsW=805 kHz) vs CISPR 32', fontsize=12)
fig.tight_layout()
fig.savefig(OUT + 'fig_emi_cispr_v2.png', dpi=140)
log('Wrote fig_emi_cispr_v2.png')

# FIG 2 : sensitivity sweep (full load, Class B)
sweeps = [
    ('10uF ESL [nH]', 'esl_10u', np.linspace(0.3e-9, 3.0e-9, 12), lambda x: x * 1e9),
    ('L1-branch trace L [nH]', 'ltr_loop', np.linspace(0, 60e-9, 13), lambda x: x * 1e9),
    ('hot-loop trace L [nH]', 'lhot', np.linspace(0, 40e-9, 13), lambda x: x * 1e9),
    ('30uF bank trace L [nH]', 'ltr_capA', np.linspace(0, 15e-9, 11), lambda x: x * 1e9),
]
fig2, axs = plt.subplots(2, 3, figsize=(16, 9))
axs = axs.ravel()
sens = {}
for k, (xl, key, vals, sc) in enumerate(sweeps):
    ys = []
    for v in vals:
        Pt = clone(Pn); Pt[key] = v
        m, *_ = worst_margin('pi', Pt, 3.0, 0.5, TR_NOM, TF_NOM)
        ys.append(m)
    sens[key] = (vals, np.array(ys))
    axs[k].plot(sc(vals), ys, '-o', color='#1f77b4')
    axs[k].set_xlabel(xl); axs[k].set_ylabel('worst Class-B margin [dB]')
    axs[k].grid(alpha=0.3); axs[k].set_title('sweep: ' + xl)
# tr sweep
trs = np.linspace(2e-9, 30e-9, 15); ys = []
for tr in trs:
    m, *_ = worst_margin('pi', Pn, 3.0, 0.5, tr, tr); ys.append(m)
sens['tr'] = (trs, np.array(ys))
axs[4].plot(trs * 1e9, ys, '-o', color='#d62728')
axs[4].set_xlabel('rise/fall time tr=tf [ns]'); axs[4].set_ylabel('worst Class-B margin [dB]')
axs[4].grid(alpha=0.3); axs[4].set_title('sweep: switching edge tr=tf')
# LISN DM impedance discrete
lisn_pts = [50.0, 100.0]; ys = []
for rr in lisn_pts:
    Pt = clone(Pn); Pt['rs_dm'] = rr
    m, *_ = worst_margin('pi', Pt, 3.0, 0.5, TR_NOM, TF_NOM); ys.append(m)
lisn_ys = ys
axs[5].bar(['50 ohm\n(v1 basis)', '100 ohm\n(strict DM)'], ys, color=['#1f77b4', '#d62728'])
axs[5].set_ylabel('worst Class-B margin [dB]'); axs[5].set_title('sweep: DM LISN impedance')
axs[5].grid(alpha=0.3, axis='y')
# log sweep results
log('')
log('--- sensitivity sweeps (full load, worst Class-B margin) ---')
for xl, key, vals, sc in sweeps:
    vv, yy = sens[key]
    log('  %-24s margin range %.1f .. %.1f dB  (over %g..%g)' %
        (xl, yy.min(), yy.max(), sc(vv).min(), sc(vv).max()))
log('  %-24s margin range %.1f .. %.1f dB  (over %g..%g ns)' %
    ('tr=tf', sens['tr'][1].min(), sens['tr'][1].max(), trs.min() * 1e9, trs.max() * 1e9))
log('  LISN DM 50ohm -> %.1f dB ; 100ohm -> %.1f dB  (delta %+.2f dB)' %
    (lisn_ys[0], lisn_ys[1], lisn_ys[1] - lisn_ys[0]))
fig2.suptitle('Sensitivity of worst-case Class-B margin (full load, nominal unless swept)', fontsize=12)
fig2.tight_layout()
fig2.savefig(OUT + 'fig_sensitivity.png', dpi=140)
log('Wrote fig_sensitivity.png')

# FIG 3 : pi vs no-pi
fig3, ax3 = plt.subplots(figsize=(11.5, 6.8))
ax3.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.2, label='CISPR 32 Class B')
ax3.plot(fgp / 1e6, cispr_A(fgp), 'k--', lw=1.3, label='CISPR 32 Class A')
for name, il, ipk, duty, col in CASES:
    f, a, z, dbv, _ = spectrum('nopi', Pnopi, ipk, duty, TR_NOM, TF_NOM)
    ax3.semilogx(f / 1e6, dbv, '--^', ms=4, color=col, alpha=0.6,
                 label='NO pi (%s, min-marg=%+.1f dB)' % (name, (cispr_B(f) - dbv).min()))
for name, il, ipk, duty, col in CASES:
    ax3.semilogx(nom[name]['f'] / 1e6, nom[name]['dbv'], '-o', ms=4, color=col,
                 label='WITH pi (%s, min-marg=%+.1f dB)' % (name, nom[name]['m'].min()))
ax3.set_xlim(0.15, 30); ax3.set_ylim(-30, 110)
ax3.set_xlabel('Frequency [MHz]'); ax3.set_ylabel('DM noise on 50 $\\Omega$ LISN [dB$\\mu$V]')
ax3.set_title('Value of the input pi-filter: DM emission WITH vs WITHOUT pi filter (3x10uF bank + L1)\n'
              'TPS62933, fsw=805 kHz, nominal parasitics')
ax3.grid(True, which='both', alpha=0.3); ax3.legend(loc='upper right', fontsize=8.5)
fig3.tight_layout()
fig3.savefig(OUT + 'fig_emi_cispr_nofilter.png', dpi=140)
log('Wrote fig_emi_cispr_nofilter.png')

# ---- numbers file ----
with open(OUT + 'emi_numbers_v2.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
log('')
log('wrote emi_numbers_v2.txt')

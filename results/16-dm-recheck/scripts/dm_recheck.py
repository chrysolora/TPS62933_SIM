#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dm_recheck.py -- robustness recheck of TPS62933 DM conducted-EMI model v2
Reproduces baseline, then: (a) geometric trace-L bounds, (b) per-param sweeps,
(c) combined worst case, (d) hidden-idealization checks (source spectrum,
QP-vs-AV limit, DC-bias derating).
Only writes to /mnt/raid10/sim-work/tps62933/dm_recheck/.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/mnt/raid10/sim-work/tps62933/dm_recheck/'
RES = '/mnt/raid10/sim-work/tps62933/results_v2/'
FSW = 805e3
np.seterr(all='ignore')

# ---------------- model functions (verbatim from emi_model_v2.py) -----------
def Zc(f, C): return 1.0 / (1j * 2 * np.pi * f * C)
def Zl(f, L): return 1j * 2 * np.pi * f * L
def capZ(f, C, esr, esl): return esr + Zl(f, esl) + Zc(f, C)

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

def solve(n, branches, inj):
    Y = np.zeros((n, n), complex); b = np.zeros(n, complex)
    for (i, j, Z) in branches:
        if abs(Z) < 1e-15: Z = 1e-12
        y = 1.0 / Z
        if i >= 0: Y[i, i] += y
        if j >= 0: Y[j, j] += y
        if i >= 0 and j >= 0: Y[i, j] -= y; Y[j, i] -= y
    b[inj] += 1.0
    return np.linalg.solve(Y, b)

def damping(P, fr):
    zc = capZ(fr, P['Ce'], P['esr_e'], P['esl_e'])
    return P['r85'] + 1.0 / (1.0 / zc + 1.0 / zc)

def Zt(net, fr, P):
    b = []
    b.append((0, -1, P['rs_dm'] + Zc(fr, P['c_lisn'])))
    b.append((0, -1, Zl(fr, P['l_mains'])))
    if net == 'pi':
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
    else:
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

def cispr_B(f):
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    66.0 - 10.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 56.0, 60.0))

def cispr_B_avg(f):
    """CISPR 32 Class B AVERAGE limit: 56->46 (150k-500k), 46 (500k-5M), 50 (5-30M)."""
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    56.0 - 10.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 46.0, 50.0))

def P_nominal():
    return dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
                esr_10u=2e-3, esl_10u=1.0e-9, esr_100n=20e-3, esl_100n=0.8e-9,
                esr_e=4.0, esl_e=3.0e-9, r85=0.1,
                L1=1e-6, dcr_l1=25e-3, cp_l1=4e-12,
                lf1=5e-9, ltr_loop=15e-9, lhot=8e-9, ltr_capA=3e-9, rf1=5e-3,
                rs_dm=50.0, c_lisn=1e-6, l_mains=50e-6)

def P_pessimistic():
    return dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
                esr_10u=6e-3, esl_10u=1.5e-9, esr_100n=40e-3, esl_100n=1.2e-9,
                esr_e=8.0, esl_e=5.0e-9, r85=0.1,
                L1=1e-6, dcr_l1=30e-3, cp_l1=6e-12,
                lf1=10e-9, ltr_loop=30e-9, lhot=20e-9, ltr_capA=8e-9, rf1=10e-3,
                rs_dm=100.0, c_lisn=1e-6, l_mains=50e-6)

def mk(p):
    P = dict(p); P['C10u'] = 10e-6; P['C100n'] = 100e-9; P['Ce'] = 10e-6
    return P

KMAX = int(30e6 / FSW)

def spectrum(net, P, ipk, duty, tr, tf, kmax=KMAX):
    f, a = trapezoid_harmonics(FSW, ipk, duty, tr, tf, kmax)
    z = Zt_vec(net, f, P)
    return f, a, z, 20 * np.log10(np.maximum(a * z, 1e-30) / 1e-6)

def worst_margin(net, P, ipk, duty, tr, tf, lim=cispr_B):
    f, a, z, dbv = spectrum(net, P, ipk, duty, tr, tf)
    m = lim(f) - dbv
    j = int(np.argmin(m))
    return m[j], f[j], dbv[j], lim(f)[j]

# ============================================================================
lines = []
def log(s=''):
    print(s); lines.append(s)

log('=' * 78)
log('DM-EMI ROBUSTNESS RECHECK  (TPS62933, CISPR 32 Class B, 150k-30MHz)')
log('=' * 78)

# ---- 0. baseline reproduce ----
Pn = mk(P_nominal()); Pp = mk(P_pessimistic())
log('')
log('[0] BASELINE REPRODUCE')
m0 = worst_margin('pi', Pn, 3.0, 0.5, 10e-9, 10e-9)
mp = worst_margin('pi', Pp, 3.0, 0.5, 20e-9, 20e-9)
log('  nominal  worst margin = %+.2f dB @ %.2f MHz (target +15.9)' % (m0[0], m0[1] / 1e6))
log('  pessim.  worst margin = %+.2f dB @ %.2f MHz (target  +0.7)' % (mp[0], mp[1] / 1e6))

# ---- 1. geometric trace inductance ----
log('')
log('[1] GEOMETRIC TRACE / LOOP INDUCTANCE  (basis)')
log('  Two standard PCB estimates (units: nH, mm):')
log('    (a) Grover strip self-L : L = 0.2*l*(ln(2l/(w+t)) + 0.5 + 0.2235*(w+t)/l)')
log('        -> no return-plane benefit; upper-ish bound for a single trace.')
log('    (b) microstrip loop over plane : L = Z0*sqrt(eeff)/c * l')
log('        -> assumes solid return plane directly beneath (2-layer: bottom GND).')
mu0 = 1.2566  # nH/mm
def grover(l, w, t=0.03):
    return 0.2 * l * (np.log(2 * l / (w + t)) + 0.5 + 0.2235 * (w + t) / l)
def microstrip(l, w, h=1.43, er=4.4):
    w_h = w / h
    eeff = (er + 1) / 2 + (er - 1) / 2 / np.sqrt(1 + 12 / w_h)
    if w_h <= 1:
        Z0 = 60 / np.sqrt(eeff) * np.log(8 * h / w + w / (4 * h))
    else:
        Z0 = 120 * np.pi / (np.sqrt(eeff) * (w / h + 1.393 + 0.667 * np.log(w / h + 1.444)))
    Lp = Z0 * np.sqrt(eeff) / 2.99792458e8 * 1e6  # nH/mm
    return Lp * l
# geometry cases for the HOT loop (VIN pin <-> local caps) and 30uF bank trace
geo = [('short/wide  (l=4, w=2)', 4, 2),
       ('typical     (l=8, w=1)', 8, 1),
       ('long/narrow (l=15, w=0.3)', 15, 0.3)]
log('  %-28s  Grover[nH]   microstrip[nH]' % 'hot-loop geometry')
for lab, l, w in geo:
    log('  %-28s  %8.1f    %8.1f' % (lab, grover(l, w), microstrip(l, w)))
L_lo, L_hi = microstrip(4, 2), grover(15, 0.3)
log('  => geometry-based L_hot plausible band  %.1f ... %.1f nH' % (L_lo, L_hi))
log('     model nominal = 8 nH (1 nH/mm x 8mm); pessimistic = 20 nH.')

# ---- 2. per-parameter sweeps (full load, D=0.5) ----
log('')
log('[2] PER-PARAMETER SWEEPS (full load Ipk=3A D=0.5, Class B QP)')
def sweep(key, vals, tr=10e-9):
    out = []
    for v in vals:
        Pt = mk(Pn); Pt[key] = v
        out.append(worst_margin('pi', Pt, 3.0, 0.5, tr, tr)[0])
    return np.array(out)

res = {}
res['lhot'] = (np.linspace(0, 60e-9, 25), None)
res['ltr_capA'] = (np.linspace(0, 40e-9, 25), None)
res['esl_10u'] = (np.linspace(0.3e-9, 3.0e-9, 15), None)
res['ltr_loop'] = (np.linspace(0, 60e-9, 13), None)
for k, (v, _) in res.items():
    res[k] = (v, sweep(k, v))
trs = np.linspace(1e-9, 30e-9, 20)
tr_y = np.array([worst_margin('pi', Pn, 3.0, 0.5, t, t)[0] for t in trs])
res['tr'] = (trs, tr_y)
lisn_v = np.array([25.0, 50.0, 75.0, 100.0])
lisn_y = sweep('rs_dm', lisn_v)
res['rs_dm'] = (lisn_v, lisn_y)
esr_v = np.linspace(1e-3, 10e-3, 10)
res['esr_10u'] = (esr_v, sweep('esr_10u', esr_v))

for k in ['lhot', 'ltr_capA', 'esl_10u', 'ltr_loop', 'tr', 'rs_dm', 'esr_10u']:
    v, y = res[k]
    if k == 'tr': sc = 1e9
    elif k in ('lhot', 'ltr_capA', 'ltr_loop', 'esl_10u'): sc = 1e9
    else: sc = 1.0
    log('  %-10s over %-16s : margin %+6.1f ... %+6.1f dB' %
        (k, '%g..%g' % (v.min() * sc, v.max() * sc), y.min(), y.max()))

# zero-margin thresholds
def zero_thresh(key, hi):
    v = np.linspace(0, hi, 200)
    y = sweep(key, v)
    idx = np.where(y < 0)[0]
    return v[idx[0]] * 1e9 if len(idx) else None
log('  zero-margin threshold:  L_hot ~ %.0f nH ;  L_tr_capA ~ %.0f nH'
    % (zero_thresh('lhot', 200e-9) or -1, zero_thresh('ltr_capA', 120e-9) or -1))

# ---- 3. combined worst case (all dominants at geometry-based worst) ----
log('')
log('[3] COMBINED WORST CASES')
combos = [
    ('nominal', dict()),
    ('pessimistic (model)', dict(lhot=20e-9, ltr_capA=8e-9, esl_10u=1.5e-9,
                                 esr_10u=6e-3, lf1=10e-9, ltr_loop=30e-9,
                                 rs_dm=100.0)),
    ('geom-worst L (lhot20,capA12)', dict(lhot=20e-9, ltr_capA=12e-9)),
    ('geom-worst L + bad ESL', dict(lhot=20e-9, ltr_capA=12e-9, esl_10u=2.0e-9)),
    ('very-worst L (lhot35,capA20)', dict(lhot=35e-9, ltr_capA=20e-9, esl_10u=2.5e-9)),
]
for lab, over in combos:
    Pt = mk(Pn); Pt.update(over)
    tr = 10e-9
    m, ft, vt, lt = worst_margin('pi', Pt, 3.0, 0.5, tr, tr)
    mav, _, vav, lav = worst_margin('pi', Pt, 3.0, 0.5, tr, tr, lim=cispr_B_avg)
    log('  %-30s QP margin %+6.1f dB @ %.2f MHz | AV margin %+6.1f dB @ %.2f MHz'
        % (lab, m, ft / 1e6, mav, ft / 1e6))

# ---- 4. hidden idealization: DC-bias derating of X5R ----------
log('')
log('[4] HIDDEN-IDEALIZATION CHECKS')
log(' (a) X5R DC-bias derating: effective cap value down to 40-70% of nominal.')
for scale in [1.0, 0.7, 0.5, 0.4]:
    Pt = mk(Pn)
    Pt['C10u'] = 10e-6 * scale  # both 10u banks scale
    m, ft, vt, lt = worst_margin('pi', Pt, 3.0, 0.5, 10e-9, 10e-9)
    log('   C_eff = %.0f%% of 10uF -> nominal QP margin %+6.1f dB @ %.2f MHz' % (scale * 100, m, ft / 1e6))

# (b) actual SPICE input source vs analytic trapezoid
log('')
log(' (b) SOURCE: analytic trapezoid vs actual SPICE i(vin) from E_full.txt')
def load_vin(fname):
    d = np.loadtxt(fname)
    t = d[:, 0]; it = d[:, 5]  # i(vin) value
    return t, it
try:
    t, iv = load_vin(RES + 'E_full.txt')
    # take last 40us (steady), resample uniform 2ns
    t0 = t[-1] - 40e-6
    sel = t >= t0
    ts = t[sel]; ivs = iv[sel]
    tu = np.arange(ts[0], ts[-1], 2e-9)
    iu = np.interp(tu, ts, ivs)
    iu = iu - iu.mean()
    n = len(iu); dt = 2e-9
    w = np.hanning(n)
    X = np.fft.rfft(iu * w); fr = np.fft.rfftfreq(n, dt)
    amp_u = 2 * np.abs(X) / w.sum()
    ks = np.arange(1, 10)
    log('   harmonic  analytic[A]  SPICE[A]   ratio[dB]')
    fa, aa = trapezoid_harmonics(FSW, 3.0, 0.5, 10e-9, 10e-9, 9)
    for k in ks:
        fi = k * FSW
        a_spice = np.interp(fi, fr, amp_u)
        r = 20 * np.log10(max(a_spice, 1e-9) / max(aa[k - 1], 1e-9))
        log('   %5d     %8.3f    %8.3f   %+7.1f' % (k, aa[k - 1], a_spice, r))
    log('   NOTE: half-load SPICE (E_half.txt) was NOT steady-state (flagged upstream).')
except Exception as e:
    log('   [source check failed: %s]' % e)

# (c) QP vs AV detector: model uses CW harmonic amplitude vs QP limit.
log('')
log(' (c) DETECTOR: model compares CW-equivalent harmonic amplitude to QP limit.')
log('     CISPR 32 Class B requires BOTH QP and AVERAGE limits (AV = QP - 10 dB).')
log('     Model does NO detector weighting; a single worst harmonic near a limit can')
log('     pass QP (+0.7 dB) yet its impulsive AV reading can be several dB lower.')
ln = np.logspace(np.log10(1.5e5), np.log10(3e7), 400)
f, a, z, dbv = spectrum('pi', mk(Pn), 3.0, 0.5, 10e-9, 10e-9)
mq = (cispr_B(f) - dbv)
ma = (cispr_B_avg(f) - dbv)
j = int(np.argmin(mq))
log('     nominal worst point %.2f MHz : QP margin %+.1f dB, AV-limit margin %+.1f dB'
    % (f[j] / 1e6, mq[j], ma[j]))

# ---- 5. figures ----
fig, axs = plt.subplots(2, 3, figsize=(16, 9)); axs = axs.ravel()
def pl(ax, key, xlab, sc=1e9, logx=False):
    v, y = res[key]
    xv = v * sc
    ax.plot(xv, y, '-o', color='#1f77b4')
    ax.axhline(0, color='r', ls='--', lw=1.2, label='0 dB (limit)')
    ax.axhspan(0, 10, color='orange', alpha=0.12)
    ax.set_xlabel(xlab); ax.set_ylabel('worst Class-B margin [dB]')
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    if logx: ax.set_xscale('log')
pl(axs[0], 'lhot', 'hot-loop trace L [nH]  (GND 20.7=model)')
axs[0].axvline(8, color='k', ls=':', label='nominal 8'); axs[0].axvline(20, color='gray', ls=':', label='pess. 20'); axs[0].legend(fontsize=7)
pl(axs[1], 'ltr_capA', '30uF bank trace L [nH]')
axs[1].axvline(3, color='k', ls=':'); axs[1].axvline(8, color='gray', ls=':')
pl(axs[2], 'esl_10u', '10uF ESL [nH]')
pl(axs[3], 'tr', 'rise/fall time tr=tf [ns]')
pl(axs[4], 'rs_dm', 'LISN DM impedance [ohm]', sc=1.0)
pl(axs[5], 'esr_10u', '10uF ESR [mohm]', sc=1e3)
fig.suptitle('DM robustness: worst Class-B margin vs each parameter (full load)', fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig_dm_robustness.png', dpi=140)
log(''); log('Wrote fig_dm_robustness.png')

with open(OUT + 'numbers_dm_recheck.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
log('wrote numbers_dm_recheck.txt')

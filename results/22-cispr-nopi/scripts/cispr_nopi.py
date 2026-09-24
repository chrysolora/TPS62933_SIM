#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cispr_nopi.py -- TPS62933 conducted CISPR 32 Class B: WITH pi vs WITHOUT pi input filter.
                 DM + CM on ONE figure (QP); AV on second panel.
Writes ONLY to /mnt/raid10/sim-work/tps62933/cispr_final_nopi/.

Caliber aligned 1:1 with cispr_final/ (same L_hot=2.70 nH real copper, same C_p=2.416 pF
real copper, same source spectrum, same CISPR 32 Class B QP+AV limits, tr=tf=10 ns).

pi-segment definition (explicit; matches board-facts.md section 2):
   pi = L1 (FXL0420-1R0-M, 1 uH) + PPHV_OUT_FILTER cap group (C73||C74 = 20 uF + C72 = 100 nF)
   The 30 uF bank (C17||C62||C63 on PPHV, BEFORE L1) is NOT part of the pi segment and is
   kept in ALL variants.
Variant (a) nopi_a : short L1 only            (both cap groups retained)
Variant (b) nopi_b : whole pi removed         (L1 shorted + C73/C74/C72 removed)

DM  model = dm_recheck/dm_recheck.py   (read-only reuse; net 'pi' verbatim)
CM  model = cm_redo/cm_redo_model.py   (read-only reuse; V_cm = 25ohm*2*pi*f*C_p*|V_sw|)
"""
import os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')

OUT = '/mnt/raid10/sim-work/tps62933/cispr_final_nopi'
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

# real copper geometry inductance (verbatim from cispr_final.py / dm_real.py)
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
L_HOT_MIC = microstrip(9.3, 3.0)      # 2.70 nH real (return-plane)
L_HOT_GRO = grover(9.3, 3.0)          # 4.44 nH real (no-return upper)
L_TR_REAL = 2.70                      # 30uF-bank segment (same geometry)

# -------- parametrized DM Zt with pi-filter modes ----------------------------
def Zt_mode(net, fr, P):
    """net in {'pi','nopi_a','nopi_b'}. Node map: 0=A(LISN), 1=P(30uF bank),
       2=V(VIN pin = source), 3=B(PPHV_OUT_FILTER cap group + damping)."""
    b = []
    b.append((0, -1, P['rs_dm'] + Zc(fr, P['c_lisn'])))
    b.append((0, -1, Zl(fr, P['l_mains'])))
    b.append((0, 1, P['rf1'] + Zl(fr, P['lf1'])))
    # 30 uF bank (C17||C62||C63) : kept in ALL variants
    for _ in range(3):
        b.append((1, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'] + P['ltr_capA'])))
    if net == 'pi':
        L1, dcr, cpl1, lloop = P['L1'], P['dcr_l1'], P['cp_l1'], P['ltr_loop']
    else:                                   # both no-pi variants: L1 shorted
        L1, dcr, cpl1, lloop = 0.0, 0.0, 0.0, 0.0
    b.append((1, 2, dcr + Zl(fr, L1 + lloop)))
    if cpl1 > 0: b.append((1, 2, Zc(fr, cpl1)))
    b.append((2, 3, Zl(fr, P['lhot'])))
    if net != 'nopi_b':                     # PPHV_OUT_FILTER group kept for pi & nopi_a
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

def worst_margin_mode(net, P, ipk, duty, tr, tf, lim=cispr_B):
    f, a, z, dbv = spectrum_mode(net, P, ipk, duty, tr, tf)
    m = lim(f) - dbv
    j = int(np.argmin(m))
    return m[j], f[j], dbv[j]

# -------- CM model functions (verbatim from cm_redo_model.py) -----------------
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

CASES = [('noload', 0.02, 0.22, 0.15, '#1f77b4'),
         ('half',   1.50, 1.50, 0.50, '#2ca02c'),
         ('full',   3.00, 3.00, 0.50, '#d62728')]
VARIANTS = [('pi', 'with-pi'), ('nopi_a', 'no-pi (a) L1 shorted'),
            ('nopi_b', 'no-pi (b) whole pi removed')]

L = []
def log(s=''):
    print(s); L.append(s)

log('='*80)
log('TPS62933 CONDUCTED CISPR 32 Class B -- WITH pi vs WITHOUT pi input filter')
log('='*80)
log('pi segment = L1 (1uH FXL0420-1R0-M) + PPHV_OUT_FILTER group (C73||C74=20uF + C72=100nF)')
log('  30uF bank (C17||C62||C63 on PPHV, BEFORE L1) is NOT in the pi segment -> kept in all variants')
log('(a) nopi_a = short L1 only  |  (b) nopi_b = L1 shorted + C73/C74/C72 removed')
log('Caliber = cispr_final: L_hot=%.2f nH real copper, C_p=%.3f pF real copper, tr=tf=10ns, CISPR32 B' % (L_HOT_MIC, CP_NOM*1e12))
log('')

# ---------------- DM per case per variant -------------------------------------
log('--- DM (real copper L_hot=%.2f nH) : 3 load cases x 3 variants ---' % L_HOT_MIC)
log('%-7s %-22s %8s %8s  %-16s %-16s' % ('case', 'variant', 'Ipk[A]', '|Zt|1', 'QP margin', 'AV margin'))
dm = {}
for nm, il, ipk, duty, col in CASES:
    for vk, vlab in VARIANTS:
        P = mk(P_nominal()); P['lhot'] = L_HOT_MIC*1e-9; P['ltr_capA'] = L_TR_REAL*1e-9
        f, a, z, dbv = spectrum_mode(vk, P, ipk, duty, 10e-9, 10e-9)
        mq = cispr_B(f) - dbv; ma = cispr_B_avg(f) - dbv
        jq = int(np.argmin(mq)); ja = int(np.argmin(ma))
        dm[(nm, vk)] = dict(f=f, dbv=dbv, mq=mq, ma=ma, col=col, il=il, jq=jq, ja=ja,
                            qp=mq[jq], av=ma[ja], fq=f[jq], fa=f[ja])
        log('%-7s %-22s %8.3f %8.1e  QP %+6.1f @%5.2fMHz  AV %+6.1f @%5.2fMHz'
            % (nm, vlab, ipk, z[0], mq[jq], f[jq]/1e6, ma[ja], f[ja]/1e6))
# conservative real copper (Grover) hot-loop L for the pi baseline, to show range
for nm, il, ipk, duty, col in CASES:
    P = mk(P_nominal()); P['lhot'] = L_HOT_GRO*1e-9; P['ltr_capA'] = L_TR_REAL*1e-9
    m, fq, _ = worst_margin_mode('pi', P, ipk, duty, 10e-9, 10e-9)
    ma_, fa_, _ = worst_margin_mode('pi', P, ipk, duty, 10e-9, 10e-9, lim=cispr_B_avg)
    log('   [pi, L_hot=%.2fnH Grover] %-6s QP %+6.1f @%5.2fMHz  AV %+6.1f dB'
        % (L_HOT_GRO, nm, m, fq/1e6, ma_))
log('')

# ---------------- CM (independent of pi in this model) ------------------------
log('--- CM (C_p=%.3f pF real copper) : 3 load cases (pi-independent in this caliber) ---' % (CP_NOM*1e12))
log('%-7s  %8s  %8s  %-18s  %-18s' % ('case', 'Vsw1[V]', 'Vcm1', 'QP margin', 'AV margin'))
cm = {}
for nm, il, ipk, duty, col in CASES:
    t, v = load_sw(nm)
    fk, am = sw_harmonics(t, v)
    d = dbuv(RCM*(2*np.pi*fk*CP_NOM*am))
    mq = cispr_B(fk) - d; ma = cispr_B_avg(fk) - d
    jq = int(np.argmin(mq)); ja = int(np.argmin(ma))
    cm[nm] = dict(f=fk, dbv=d, mq=mq, ma=ma, col=col, il=il, jq=jq, ja=ja,
                  qp=mq[jq], av=ma[ja], fq=fk[jq], fa=fk[ja])
    log('%-7s  %8.3f  %8.1f  QP %+6.1f @%5.2fMHz  AV %+6.1f @%5.2fMHz %s'
        % (nm, am[0], d[0], mq[jq], fk[jq]/1e6, ma[ja], fk[ja]/1e6,
           '' if nm != 'half' else '(half NOT steady-state)'))
log('  NOTE: pi-segment elements (L1 + caps) all return to on-board GND; the CM loop is')
log('        SW -> C_p(board<->earth) -> LISN(25ohm) -> cable, so no pi element lies in it')
log('        -> CM curve is IDENTICAL for pi / nopi_a / nopi_b in this caliber.')
log('')

# ---------------- verification: with-pi reproduces cispr_final ----------------
log('--- VERIFY with-pi reproduces cispr_final/ (DM full +25.1/+15.1, CM full -16.0) ---')
log('  DM  full  QP %+.1f  AV %+.1f   (target +25.1 / +15.1)' % (dm[('full','pi')]['qp'], dm[('full','pi')]['av']))
log('  CM  full  QP %+.1f  AV %+.1f   (target -16.0 / -26.0)' % (cm['full']['qp'], cm['full']['av']))
log('')

with open(OUT + '/numbers_nopi.txt', 'w') as fh:
    fh.write('\n'.join(L) + '\n')

# ---------------- FIGURE: QP (DM+CM) + AV -------------------------------------
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)
fig, (axQ, axA) = plt.subplots(1, 2, figsize=(17.5, 7.4))
lstyle = {'pi': '-', 'nopi_a': '--', 'nopi_b': '-.'}
for ax, lim, ttl in ((axQ, cispr_B, 'QUASI-PEAK'), (axA, cispr_B_avg, 'AVERAGE')):
    ax.plot(fgp/1e6, lim(fgp), 'k-', lw=2.6, label='CISPR 32 Class B %s limit' % ttl)
    ax.set_xscale('log'); ax.set_xlim(0.15, 30); ax.set_ylim(0, 110)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel('LISN voltage [dB$\\mu$V]')
    ax.set_title('%s detector' % ttl)
# DM curves: colour by load case, linestyle by variant
for nm, il, ipk, duty, col in CASES:
    for vk, vlab in VARIANTS:
        r = dm[(nm, vk)]
        axQ.semilogx(r['f']/1e6, r['dbv'], lstyle[vk], color=col, lw=1.5, alpha=0.9,
                     label='DM %-6s %-17s QPm %+.1f' % (nm, vlab, r['qp']))
        axA.semilogx(r['f']/1e6, r['dbv'], lstyle[vk], color=col, lw=1.5, alpha=0.9,
                     label='DM %-6s %-17s AVm %+.1f' % (nm, vlab, r['av']))
# CM curves: markers, one per case (identical for all pi variants)
for nm, il, ipk, duty, col in CASES:
    r = cm[nm]
    for j in range(0, len(r['f']), 3):
        pass
    axQ.semilogx(r['f']/1e6, r['dbv'], ':', color=col, lw=2.2, marker='o', ms=2.5,
                 label='CM %-6s Ild=%.1fA  QPm %+.1f' % (nm, il, r['qp']))
    axA.semilogx(r['f']/1e6, r['dbv'], ':', color=col, lw=2.2, marker='o', ms=2.5,
                 label='CM %-6s Ild=%.1fA  AVm %+.1f' % (nm, il, r['av']))
axQ.text(0.16, 105, 'DM: solid=with-pi  dashed=no-pi(a) L1short  dashdot=no-pi(b) pi-removed\n'
                    'CM: dotted+o (identical for all pi variants in this caliber)',
         fontsize=8.0, va='top')
axQ.legend(loc='lower left', fontsize=6.6, ncol=2)
axA.legend(loc='lower left', fontsize=6.6, ncol=2)
fig.suptitle('TPS62933 conducted emission vs CISPR 32 Class B -- WITH vs WITHOUT pi input filter\n'
             'fsw=805 kHz | DM: L_hot=2.70 nH real copper | CM: C_p=2.416 pF real copper | caliber = cispr_final/',
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT + '/fig_cispr_nopi_compare.png', dpi=150)
log('Wrote fig_cispr_nopi_compare.png')

with open(OUT + '/numbers_nopi.txt', 'w') as fh:
    fh.write('\n'.join(L) + '\n')

# ---------------- comparison table (markdown) ---------------------------------
def g(x):   # grade
    return 'B'
def cell(v, f_):
    return '**%+.1f** @%.2f' % (v, f_/1e6)
rows = []
rows.append('# comparison_table.md — TPS62933 CISPR 32 Class B 裕量对比（有π / 无π-a / 无π-b）\n')
rows.append('**判据**：裕量 = 限值 − 发射。**正=过**，**负=超**。每格标 A/B/C。\n')
rows.append('**π 段定义**：`L1(1µH) + PPHV_OUT_FILTER 电容组(C73‖C74=20µF + C72=100nF)`；')
rows.append('30µF 银行（PPHV 上、L1 前）不属 π 段，**三档全部保留**。\n')
rows.append('- 无π-a = 只短路 L1（两组电容全保留）\n- 无π-b = 整个 π 段移除（L1 短路 + 去掉 C73/C74/C72）\n')
rows.append('口径与 `cispr_final/` 一致：L_hot=2.70 nH（真实铜皮）、C_p=2.416 pF（真实铜皮）、')
rows.append('源谱 results_v2/E_*.txt v(sw) 与解析梯形波(tr=tf=10ns)、CISPR 32 Class B QP+AV。\n')

rows.append('## DM（差模，L_hot=2.70 nH）\n')
rows.append('| 工况 | 有π QP | 有π AV | 无π-a QP | 无π-a AV | 无π-b QP | 无π-b AV | 等级 |')
rows.append('|---|---|---|---|---|---|---|---|')
for nm, il, ipk, duty, col in CASES:
    a = dm[(nm,'pi')]; b = dm[(nm,'nopi_a')]; c = dm[(nm,'nopi_b')]
    rows.append('| %s (%.2f A) | %s | %s | %s | %s | %s | %s | B |'
                % (nm, il, cell(a['qp'],a['fq']), cell(a['av'],a['fa']),
                   cell(b['qp'],b['fq']), cell(b['av'],b['fa']),
                   cell(c['qp'],c['fq']), cell(c['av'],c['fa'])))

rows.append('\n## CM（共模，C_p=2.416 pF）— 本口径下与 π 无关，三档相同\n')
rows.append('| 工况 | 有π QP | 有π AV | 无π-a QP | 无π-a AV | 无π-b QP | 无π-b AV | 等级 |')
rows.append('|---|---|---|---|---|---|---|---|')
for nm, il, ipk, duty, col in CASES:
    r = cm[nm]; tag = 'C*' if nm == 'half' else 'B'
    same = '（同）'
    rows.append('| %s (%.2f A) | %s | %s | %s | %s | %s | %s | %s |'
                % (nm, il, cell(r['qp'],r['fq']), cell(r['av'],r['fa']),
                   '（同）', '（同）', '（同）', '（同）', tag))

rows.append('\n## π 滤波器增益（DM 满载，限值−发射的差值 = 有π − 无π）\n')
rows.append('| 对比 | QP Δ | AV Δ |')
rows.append('|---|---|---|')
af = dm[('full','pi')]; bf = dm[('full','nopi_a')]; cf = dm[('full','nopi_b')]
rows.append('| 有π − 无π-a | %+.1f dB | %+.1f dB |' % (af['qp']-bf['qp'], af['av']-bf['av']))
rows.append('| 有π − 无π-b | %+.1f dB | %+.1f dB |' % (af['qp']-cf['qp'], af['av']-cf['av']))
rows.append('\n（正值 = π 段带来的裕量提升；＝ 无π发射 − 有π发射，同为 Zt 比值 dB。）\n')

rows.append('## 等级说明\n')
rows.append('- **B**：趋势/相对可信，**绝对值不认证**——DM 绝对值依赖假定 SW 边沿 tr=tf=10 ns（解析梯形波）')
rows.append('  与电路仿真来源的 ipk/duty；CM 依赖集总 V_cm=25Ω·2πf·C_p·|V_sw| 近似。与 `cispr_final/` 分级一致。')
rows.append('- **C\\***：半载 SPICE 未达稳态（沿用上游标记），仅作参考。')
rows.append('- 关键量 C_p=2.416 pF、L_hot=2.70 nH 本身为 **A**（`cispr_final/` 已交叉验证）。')

with open(OUT + '/comparison_table.md', 'w') as fh:
    fh.write('\n'.join(rows) + '\n')
log('Wrote comparison_table.md')
print('DONE')

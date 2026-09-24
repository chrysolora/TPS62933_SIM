# -*- coding: utf-8 -*-
# cm_fix: 传导 CM -- 修「耦合元件用错 + 回路分配未定」
# Model A: 一般做法 C_p = 开关节点->大地寄生, 全部经 LISN (无 f_exit)
# Model B: cm_audit 的 SW<->板内GND 电容 C_sw=2.416pF, 仅 f_exit 部分外流, f_exit 用阻抗分压真算
import os, math
import numpy as np

BASE = "/mnt/raid10/sim-work/tps62933"
OUT = BASE + "/cm_fix"
os.makedirs(OUT, exist_ok=True)
_lines = []
def L(*a):
    s = " ".join(str(x) for x in a); _lines.append(s); print(s, flush=True)

FSW = 805e3; M = 32; KMAX = int(30e6 / FSW)
A_BOARD = 0.0265 * 0.0505
CP_SWG = 2.416e-12
CY = 2.2e-9
L_CAB = 100e-9
R_LISN = 25.0
ZLOC = [1.0, 2.0, 6.0]

def load_sw(case):
    a = np.loadtxt("%s/results_v2/E_%s.txt" % (BASE, case))
    return a[:, 0], a[:, 3]

def sw_harmonics(t, vsw, M=M):
    Tend = t[-1]; Tw = M / FSW; t0 = Tend - Tw
    m = t >= t0; tt = t[m]; vv = vsw[m]
    nu = int(round(Tw * FSW * 128))
    tu = np.linspace(t0, Tend, nu); vv_u = np.interp(tu, tt, vv)
    V = np.fft.rfft(vv_u - vv_u.mean()) / len(vv_u) * 2.0
    df = 1.0 / Tw
    ks = np.arange(1, KMAX + 1); tg = ks * FSW
    idx = np.clip(np.round(tg / df).astype(int), 0, len(V) - 1)
    return tg, np.abs(V[idx])

def z_out(f):
    w = 2 * np.pi * np.maximum(f, 1.0)
    zy = 1.0 / (1j * w * CY)
    zc = 1j * w * L_CAB
    return 1.0 / (1.0 / R_LISN + 1.0 / (zc + zy))

def f_exit(f, Zloc):
    zo = np.abs(z_out(f))
    return Zloc / (Zloc + zo)

def lim_qp(f):
    fm = np.asarray(f, dtype=float) / 1e6
    return np.where(fm < 0.5, 66.0 - 19.1 * np.log10(np.maximum(fm, 1e-9) / 0.15),
           np.where(fm < 5.0, 56.0, 60.0))

def margin(Vcm_dBuV, f):
    return lim_qp(f) - Vcm_dBuV

def vu(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-30) / 1e-6)

sw = {}
L("=" * 74)
L("TPS62933 CM-fix : corrected coupling element + REAL f_exit")
L("=" * 74)
a = np.loadtxt("%s/results_v2/E_full.txt" % BASE)
L("col sanity (full, amplitude over last 32 periods; t_end=%.1f us):" % (a[-1, 0] * 1e6))
for c in range(1, a.shape[1]):
    tg, am = sw_harmonics(a[:, 0], a[:, c])
    L("  col %2d : |Vsw1|=%.3f  |Vsw3|=%.4f  |Vsw5|=%.4f" % (c, am[0], am[2], am[4]))
L("-> using col 3 as v(sw) (consistent with cm_redo)")

for case in ["noload", "half", "full"]:
    t, v = load_sw(case); tg, am = sw_harmonics(t, v)
    sw[case] = (tg, am)
    L("case %-7s |Vsw1|=%.3f  |Vsw5|=%.3f  |Vsw37|=%.3f" % (case, am[0], am[4], am[36]))
L("")

L("--- C_be (board<->earth parallel-plate, A=26.5x50.5mm) ---")
for d in [0.01, 0.10, 0.40]:
    L("  d=%.0f cm : C_be=%.4f pF" % (d * 100, 8.854e-12 * A_BOARD / d * 1e12))
L("")

L("=== Model A: general practice (C_p = SW-to-earth stray, ALL through LISN) ===")
L("  case    C_p[pF]  worst-QP[dB] @f        worst-AV[dB] @f")
modelA = {}
for case in ["noload", "half", "full"]:
    tg, am = sw[case]; modelA[case] = {}
    for cp in [1.0, 2.416, 5.0, 10.0, 20.0]:
        Vcm = R_LISN * 2 * np.pi * tg * (cp * 1e-12) * am
        dcm = vu(Vcm); mq = margin(dcm, tg); k = int(np.argmin(mq))
        ma = mq - 10.0; k2 = int(np.argmin(ma))
        modelA[case][cp] = (mq[k], tg[k] / 1e6, ma[k2], tg[k2] / 1e6)
        L("  %-7s %7.3f   %8.1f @%8.3fM   %8.1f @%8.3fM" % (case, cp, mq[k], tg[k] / 1e6, ma[k2], tg[k2] / 1e6))
L("")

L("=== Model B: cm_audit SW<->board-GND C_sw=2.416pF, REAL exit fraction f_exit ===")
L("  f_exit = |Z_loc| / (|Z_loc| + |Z_out|), Z_out = LISN(25ohm) || (cable + Ycap 2.2nF)")
modelB = {}
for case in ["noload", "half", "full"]:
    tg, am = sw[case]; modelB[case] = {}
    for zl in ZLOC:
        fe = f_exit(tg, zl)
        Vcm = R_LISN * fe * 2 * np.pi * tg * CP_SWG * am
        dcm = vu(Vcm); mq = margin(dcm, tg); k = int(np.argmin(mq))
        ma = mq - 10.0; k2 = int(np.argmin(ma))
        modelB[case][zl] = (mq[k], tg[k] / 1e6, ma[k2], tg[k2] / 1e6)
        fe1 = 100 * fe[np.argmin(np.abs(tg - FSW))]
        L("  %-7s Zloc=%4.1fohm  f_exit@0.805M=%.1f%%  worst-QP=%7.1f @%.3fM  worst-AV=%7.1f @%.3fM"
          % (case, zl, fe1, mq[k], tg[k] / 1e6, ma[k2], tg[k2] / 1e6))
L("")
L("f_exit(0.805MHz): " + ", ".join("Zloc=%.1f->%.1f%%" % (z, 100 * f_exit(np.array([FSW]), z)[0]) for z in ZLOC))
L("f_exit(30MHz)   : " + ", ".join("Zloc=%.1f->%.1f%%" % (z, 100 * f_exit(np.array([30e6]), z)[0]) for z in ZLOC))
L("")

allA = [modelA[c][cp][0] for c in modelA for cp in modelA[c]]
allB = [modelB[c][z][0] for c in modelB for z in modelB[c]]
L("=== VERDICT ===")
L(" Model A QP-margin range: %.1f .. %.1f dB -> %s" % (min(allA), max(allA), "ALL EXCEED" if max(allA) < 0 else "some pass"))
L(" Model B QP-margin range: %.1f .. %.1f dB -> %s" % (min(allB), max(allB), "ALL EXCEED" if max(allB) < 0 else "some/most PASS"))

try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    tg, am = sw["full"]; fm = tg / 1e6
    A = R_LISN * 2 * np.pi * tg * 2.416e-12 * am
    ax[0].semilogx(fm, vu(A), "r-", label="Model A: C_p=2.416pF (all exits)")
    for zl, sty in zip(ZLOC, ["b-", "g-", "k-"]):
        fe = f_exit(tg, zl); B = R_LISN * fe * 2 * np.pi * tg * CP_SWG * am
        ax[0].semilogx(fm, vu(B), sty, label="Model B: Csw=2.416pF f_exit(Zloc=%.0f)" % zl)
    lq = lim_qp(tg)
    ax[0].semilogx(fm, lq, "m--", label="CISPR32 B QP")
    ax[0].semilogx(fm, lq - 10, "c:", label="CISPR32 B AV")
    ax[0].set_ylabel("V_cm (dBuV)"); ax[0].set_ylim(0, 120)
    ax[0].legend(fontsize=7); ax[0].grid(True, which="both", alpha=.3)
    ax[0].set_title("TPS62933 conducted CM (full load) - Model A vs Model B")
    ax[1].semilogx(fm, 100 * f_exit(tg, 2.0), "b-")
    ax[1].set_ylabel("f_exit @Zloc=2ohm (%)"); ax[1].set_xlabel("f (MHz)")
    ax[1].grid(True, which="both", alpha=.3)
    plt.tight_layout(); plt.savefig(OUT + "/fig_cm_fix.png", dpi=130)
    L("plot: " + OUT + "/fig_cm_fix.png")
except Exception as e:
    L("plot FAILED: %r" % e)

open(OUT + "/numbers_cm_fix.txt", "w").write("\n".join(_lines) + "\n")
print("WROTE", OUT + "/numbers_cm_fix.txt")

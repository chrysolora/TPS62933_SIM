# -*- coding: utf-8 -*-
# dm_real: 用本板【真实铺铜几何】定 DM 的两段电感，重算 CISPR 32 QP/AV
# 复用 dm_recheck.py 的模型(dm_recheck.py 顶层 prefix 不过 exec 到 lines=[])
import numpy as np, sys
SRC = "/mnt/raid10/sim-work/tps62933/dm_recheck/dm_recheck.py"
src = open(SRC).read()
prefix = src.split("lines = []")[0]
G = {}
exec(prefix, G)
P_nominal = G["P_nominal"]; mk = G["mk"]
worst_margin = G["worst_margin"]; cispr_B = G["cispr_B"]; cispr_B_avg = G["cispr_B_avg"]

def grover(l, w, t=0.03):
    return 0.2 * l * (np.log(2 * l / (w + t)) + 0.5 + 0.2235 * (w + t) / l)
def microstrip(l, w, h=1.43, er=4.4):
    w_h = w / h
    eeff = (er + 1) / 2 + (er - 1) / 2 / np.sqrt(1 + 12 / w_h)
    if w_h <= 1:
        Z0 = 60 / np.sqrt(eeff) * np.log(8 * h / w + w / (4 * h))
    else:
        Z0 = 120 * np.pi / (np.sqrt(eeff) * (w / h + 1.393 + 0.667 * np.log(w / h + 1.444)))
    return Z0 * np.sqrt(eeff) / 2.99792458e8 * 1e6 * l

OUT = "/mnt/raid10/sim-work/tps62933/dm_real"
import os; os.makedirs(OUT, exist_ok=True)
_log = []
def L(*a):
    s = " ".join(str(x) for x in a); _log.append(s); print(s, flush=True)

L("="*76)
L("DM with THIS BOARD's REAL COPPER GEOMETRY  (TPS62933, CISPR 32 Class B)")
L("="*76)
# real copper measured from epru (fea3 parse): 24VIN top pour POUR1 = 3.0 x 9.3 mm, area 27.8mm2
L("Real copper (from pourSim.epru):")
L("  POUR1 24VIN (top): bbox 3.0 x 9.3 mm , area 27.8 mm^2  (wide POUR, not a thin trace)")
L("  POUR4 GND  (top): 31.1 x 52.8 mm  /  POUR6 GND (bot): 30.0 x 54.9 mm  (solid planes)")
L("  -> the VIN hot loop and the 30uF bank sit on the same wide 24VIN pour")
L("     whose return is the bottom GND plane at 1.43 mm.")
L("")
L("Inductance of the REAL 24VIN seg (l=9.3mm, w=3.0mm):")
gro = grover(9.3, 3.0); mic = microstrip(9.3, 3.0)
L("   Grover strip (no return benefit, upper)  = %.2f nH" % gro)
L("   microstrip over plane (w=3,h=1.43,er4.4) = %.2f nH" % mic)
L("   => real hot-loop L ~ %.1f - %.1f nH   (model nominal used 8 nH ; pess 20 nH)"
  % (min(gro, mic), max(gro, mic)))
L("   => real 30uF-bank seg similar -> ~ %.1f - %.1f nH  (model nominal 3 nH)"
  % (min(microstrip(4,3), grover(4,3)), max(microstrip(9.3,3), grover(9.3,3))))
L("")

# sweep lhot with the REAL geometry value highlighted
scan = np.linspace(0, 25e-9, 26)
qp = []; av = []
for v in scan:
    P = mk(P_nominal()); P["lhot"] = v
    qp.append(worst_margin("pi", P, 3.0, 0.5, 10e-9, 10e-9, lim=cispr_B)[0])
    av.append(worst_margin("pi", P, 3.0, 0.5, 10e-9, 10e-9, lim=cispr_B_avg)[0])
qp = np.array(qp); av = np.array(av)

def report(tag, lhot, ltr):
    P = mk(P_nominal()); P["lhot"] = lhot; P["ltr_capA"] = ltr
    mq, fq, dbq, lq = worst_margin("pi", P, 3.0, 0.5, 10e-9, 10e-9, lim=cispr_B)
    ma, fa, dba, la = worst_margin("pi", P, 3.0, 0.5, 10e-9, 10e-9, lim=cispr_B_avg)
    L("  %-34s QP %+6.1f dB @%5.2fMHz | AV %+6.1f dB @%5.2fMHz"
      % (tag, mq, fq/1e6, ma, fa/1e6))

L("=== margins with REAL geometry (vs CISPR32 Class B) ===")
report("REAL copper: lhot=%.1f,ltr=2.7nH" % mic, mic*1e-9, 2.7e-9)
report("REAL copper: lhot=%.1f,ltr=2.7nH" % ((gro+mic)/2), (gro+mic)/2*1e-9, 2.7e-9)
report("REAL copper: lhot=%.1f (conserv.)" % gro, gro*1e-9, 2.7e-9)
report("[ref] model nominal 8/3 nH", 8e-9, 3e-9)
report("[ref] geom worst 15/8 nH", 15e-9, 8e-9)
report("[ref] model pessimistic 20/8", 20e-9, 8e-9)
L("")
L("=== verdict (real copper) ===")
P = mk(P_nominal()); P["lhot"] = mic*1e-9; P["ltr_capA"] = 2.7e-9
mq = worst_margin("pi", P, 3.0, 0.5, 10e-9, 10e-9, lim=cispr_B)[0]
ma = worst_margin("pi", P, 3.0, 0.5, 10e-9, 10e-9, lim=cispr_B_avg)[0]
L(" REAL copper (lhot=%.1fnH): QP %+.1f dB, AV %+.1f dB -> %s"
  % (mic, mq, ma, "PASS both" if min(mq, ma) > 0 else ("QP pass, AV fail" if mq > 0 else "FAIL")))
L(" previous 'pessimistic envelope' used lhot=20nH = %.1fx the real-copper value -> that is the 'bad layout' overshoot."
  % (20.0/mic))

try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(scan*1e9, qp, "b-", label="QP margin (Class B QP)")
    ax.plot(scan*1e9, av, "r-", label="AV margin (Class B AV = QP-10dB)")
    ax.axhline(0, color="k", lw=.8)
    ax.axvline(mic, color="g", ls="--", label="REAL copper = %.1f nH" % mic)
    ax.axvline(8, color="k", ls=":", label="model nominal 8")
    ax.axvline(20, color="gray", ls=":", label="model pessimistic 20")
    ax.set_xlabel("hot-loop inductance lhot [nH]"); ax.set_ylabel("margin to limit [dB]")
    ax.set_title("DM margin vs hot-loop L - real copper is far left (better)")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    plt.tight_layout(); plt.savefig(OUT+"/fig_dm_real.png", dpi=130)
    L("plot: "+OUT+"/fig_dm_real.png")
except Exception as e:
    L("plot failed: %r" % e)

open(OUT+"/numbers_dm_real.txt","w").write("\n".join(_log)+"\n")
print("WROTE", OUT+"/numbers_dm_real.txt")

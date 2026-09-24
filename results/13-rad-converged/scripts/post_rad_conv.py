#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""post_rad_conv.py -- converged radiation model for the TPS62933 copper cutout question.

Fixes the B2 non-convergence by (a) removing geometry-induced micro-cells (all
geometry quantized to a 0.25 mm FDTD grid -> dt ~12x larger), (b) running long
enough that the Gaussian excitation completes AND the field decays (EndCriteria
1e-4 = -40 dB), and (c) keeping the physical normalisation  H = E_norm / |if_tot|
(the lumped port current IS the physical hot-loop current, in series with the loop).

Reads ver_<v>_dm_<L>/ {nf2ff.npz,port.npz}. Writes only to rad_conv/.
"""
import os, re, numpy as np, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
VERS  = ["full", "top", "dual"]
LABEL = {"full": "FullCopper", "top": "TopCutout", "dual": "DualCutout"}
COL   = {"full": "C0", "top": "C2", "dual": "C3"}
LEVELS = [("C0", 1.5), ("C1", 1.0), ("C2", 0.75), ("C3", 0.5)]     # coarse -> fine
DM_FSW, DM_D, DM_TR = 805e3, 0.5, 5e-9
IPK = 4.3                                              # 3.0 A load
FT = [30e6, 100e6, 300e6]

def harm(fmin, fmax, fsw, A, D, tr):
    N = 400001
    t = np.linspace(0, 1.0/fsw, N, endpoint=False)
    T = 1.0/fsw; top = max(D*T - 2*tr, 0.0)
    x = np.interp(t, [0, tr, tr+top, 2*tr+top, T], [0, A, A, 0, 0])
    X = np.abs(np.fft.rfft(x)/N)
    k = np.arange(1, int(fmax/fsw)+1); f = k*fsw
    m = (f >= fmin) & (f <= fmax)
    return f[m], X[k[m]]

def interp_cplx(H, fsrc, fdst):
    nf, a, b = H.shape; Hf = H.reshape(nf, -1)
    out = np.empty((len(fdst), a*b), complex)
    for c in range(a*b):
        out[:, c] = np.interp(fdst, fsrc, Hf[:, c].real) + 1j*np.interp(fdst, fsrc, Hf[:, c].imag)
    return out.reshape(len(fdst), a, b)

def dbmax(E):
    return 20*np.log10(np.maximum(np.abs(E).max(axis=(1, 2)), 1e-30)/1e-6)

def band(dB, FAP, t, half=2):
    i = int(np.argmin(np.abs(FAP - t)))
    return dB[max(0, i-half):i+half+1].max()

fh_dm, Ih_unit = harm(30e6, 1e9, DM_FSW, 1.0, DM_D, DM_TR)
FAP = fh_dm; IH = Ih_unit*IPK

def load(v, lv):
    p = os.path.join(HERE, "ver_%s_dm_%s" % (v, lv))
    if not os.path.exists(os.path.join(p, "nf2ff.npz")): return None
    d = np.load(os.path.join(p, "nf2ff.npz")); pd = np.load(os.path.join(p, "port.npz"))
    ift = np.abs(pd["dm_if_tot"])
    H = interp_cplx(d["E_norm"]/ift[:, None, None], d["freq"], FAP)
    Edm = H * IH[:, None, None] / 3.0
    return dict(dB=dbmax(Edm), E=Edm, ift=ift)

def energy_db(logname):
    try:
        for ln in reversed(open(os.path.join(HERE, logname), errors="ignore").read().splitlines()):
            if "Energy:" in ln and "dB)" in ln:
                return float(ln.split("(")[-1].split("dB")[0].replace("-", "-").replace(" ", ""))
    except Exception: return None
    return None

def conv_flag(logname):
    try:
        t = open(os.path.join(HERE, logname), errors="ignore").read()
        return ("end-criteria" in t and "of -40dB was reached" in t and "was reached before" not in t)
    except Exception: return None

lines = []
def P(s=""):
    print(s); lines.append(s)

P("=== TPS62933 rad_conv : converged board radiation (DM, Ipk=4.3A @3.0A) ===")
P("normalisation H = E_norm/|if_tot| (physical hot-loop current); E_3m = H*|I(f)|/3")
P("geometry quantized to 0.25mm grid (kills micro-cells); FR4 1.50mm; limit CISPR32-B@3m")
P("")

# ---------- convergence table ----------
P("--- E@3m [dBuV/m] vs mesh level (band-max +/-2 harmonics) ---")
data = {}
for v in VERS:
    data[v] = {}
    for lv, resb in LEVELS:
        r = load(v, lv)
        if r is None: continue
        data[v][lv] = r
        row = [band(r["dB"], FAP, t) for t in FT]
        en = energy_db("C%s_%s.log" % (lv[0], v)) if False else None
        P("  %-5s %s (resb %.2f) : E30=%6.2f  E100=%6.2f  E300=%6.2f" % (v, lv, resb, row[0], row[1], row[2]))

# end-state energy from logs
P("")
P("--- openEMS end-state energy & convergence flag (from run logs) ---")
for v in VERS:
    for lv, _ in LEVELS:
        fn = "C%s_%s.log" % (lv[0], v)
        p = os.path.join(HERE, fn)
        if not os.path.exists(p): continue
        txt = open(p, errors="ignore").read()
        m = re.findall(r"Energy:~?\s*[0-9.eE+\-]+\s*\(\s*([-+0-9.]+)\s*dB\)", txt)
        e = float(m[-1]) if m else None
        ok = ("of -40.00dB was reached" in txt and "was reached before" not in txt)
        P("  %-6s %s : final energy %s dB  end-criteria(-40dB) met: %s" %
          (fn, "", ("%+.2f" % e) if e is not None else "?", ok))

# ---------- convergence deltas ----------
P("")
P("--- convergence: |E30(level) - E30(next)| (fine pair is the acceptance test) ---")
conv_err = {}
for v in VERS:
    lvs = [lv for lv, _ in LEVELS if lv in data[v]]
    vals = [band(data[v][lv]["dB"], FAP, 30e6) for lv in lvs]
    if len(vals) >= 2:
        P("  %-5s E30 levels = %s" % (v, " ".join("%.2f" % x for x in vals)))
        for i in range(1, len(vals)):
            P("      |dE30| %s->%s = %.2f dB" % (lvs[i-1], lvs[i], abs(vals[i]-vals[i-1])))
        conv_err[v] = max(abs(vals[i]-vals[i-1]) for i in range(1, len(vals)))
    else:
        P("  %-5s only one level available" % v)

# error band = max convergence drift over the finest available pair, per variant
P("")
P("--- resolvability: three-version spread vs mesh error band ---")
# pick finest common level
common = [lv for lv, _ in LEVELS if all(lv in data[v] for v in VERS)]
if not common:
    common = [lv for lv, _ in LEVELS if any(lv in data[v] for v in VERS)]
lvc = common[-1]
# per-variant error band from its own level drift (fallback: global max)
global_band = max(conv_err.values()) if conv_err else 2.0
P("  comparison level = %s ; global convergence band = %.2f dB" % (lvc, global_band))
resolv = {}
for t in FT:
    e30 = {v: band(data[v][lvc]["dB"], FAP, t) for v in VERS}
    spread = max(e30.values()) - min(e30.values())
    band_v = max(conv_err.get(v, global_band) for v in VERS)
    P("  @%.0f MHz : Full=%.2f Top=%.2f Dual=%.2f  spread(Dual-Full)=%+.2f  band=%.2f -> %s" %
      (t/1e6, e30["full"], e30["top"], e30["dual"], e30["dual"]-e30["full"], band_v,
       "RESOLVABLE" if abs(e30["dual"]-e30["full"]) > band_v else "not resolvable (<noise)"))
    resolv[t] = (spread, band_v, e30["dual"]-e30["full"])

open(os.path.join(HERE, "numbers_conv.txt"), "w").write("\n".join(lines)+"\n")

# ---------- figures ----------
LVX = {"C0": 0, "C1": 1, "C2": 2, "C3": 3}
fig, ax = plt.subplots(figsize=(8.5, 5.4))
for v in VERS:
    lvs = [lv for lv, _ in LEVELS if lv in data[v]]
    xs = [LVX[lv] for lv in lvs]
    ys = [band(data[v][lv]["dB"], FAP, 30e6) for lv in lvs]
    ax.plot(xs, ys, "o-", color=COL[v], label=LABEL[v])
    if len(ys) >= 2:
        ax.annotate("%.1f dB" % abs(ys[-1]-ys[-2]), (xs[-1], ys[-1]),
                    textcoords="offset points", xytext=(-4, 8), fontsize=8, color=COL[v])
ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["C0\n1.5mm", "C1\n1.0mm", "C2\n0.75mm", "C3\n0.5mm"])
ax.set_xlabel("mesh level (base cell size)"); ax.set_ylabel("E @3m @30 MHz [dBuV/m]")
ax.set_title("Mesh convergence of E@3m@30MHz (DM, 3.0A)")
ax.grid(alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_rad_conv.png"), dpi=130)
print("saved fig_rad_conv.png")

fig2, ax = plt.subplots(figsize=(10, 5.6))
LIM = lambda f: np.where(f < 230e6, 40.0, 47.0)
for v in VERS:
    if lvc in data[v]:
        ax.semilogx(FAP/1e6, data[v][lvc]["dB"], color=COL[v], lw=1.0, label="%s (%s)" % (LABEL[v], lvc))
ax.semilogx(FAP/1e6, LIM(FAP), "k--", lw=1.6, label="CISPR 32 B @3m")
ax.set_xlabel("Freq [MHz]"); ax.set_ylabel("E @3m [dBuV/m]")
ax.set_title("Converged radiated far field, three copper-cutout variants (DM, 3.0A)")
ax.set_ylim(0, 60); ax.grid(True, which="both", alpha=0.3); ax.legend(loc="lower left", fontsize=8)
fig2.tight_layout(); fig2.savefig(os.path.join(HERE, "fig_rad_versions_final.png"), dpi=130)
print("saved fig_rad_versions_final.png")
print("done")

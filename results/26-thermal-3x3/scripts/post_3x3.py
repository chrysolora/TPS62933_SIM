#!/usr/bin/env python3
"""TPS62933 thermal 3x3 post-processing: figures + comparison json.
Usage: post_3x3.py
"""
import os, json, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = "/mnt/raid10/sim-work/tps62933/thermal_3x3"
os.chdir(ROOT)
COMPS = {"U21": (2.10, 1.60, 8.25, 26.88), "D1": (4.30, 2.60, 14.10, 6.56),
         "L2": (6.50, 6.50, 6.22, 33.74), "L1": (4.20, 4.20, 16.38, 18.75)}
OX, OY = 6.223, 33.736
KEEPOUT = [("L2", OX-3.2999934, OY-3.499993, OX+3.2999934, OY+3.499993),
           ("L1", OX+7.9838042, OY-17.272, OX+12.3361958, OY-12.7098044)]
GEOS = ["fullcu", "topcut", "dualcut"]
GLAB = {"fullcu": "FullCopper", "topcut": "TopCutout", "dualcut": "DualCutout"}
LOADS = [3.0, 1.5, 0.6]
LKEY = {3.0: "L30", 1.5: "L15", 0.6: "L06"}


def read_vtu(casedir):
    fs = glob.glob(os.path.join(casedir, "mesh", "*.vtu")) or glob.glob(os.path.join(casedir, "*.vtu"))
    def num(f):
        m = re.search(r"_t(\d+)\.vtu$", f); return int(m.group(1)) if m else -1
    f = sorted(fs, key=num, reverse=True)[0]
    raw = open(f, "rb").read(); txt = raw.decode("latin1")
    npts = int(re.search(r'NumberOfPoints="(\d+)"', txt).group(1))
    base = txt.index("_", txt.index("<AppendedData")) + 1
    def arr(off, dt, n):
        p = base + off
        ln = int(np.frombuffer(raw[p:p+4], dtype="<u4")[0])
        return np.frombuffer(raw[p+4:p+4+ln], dtype=dt)[:n]
    toff = int(re.search(r'Name="temperature"[^>]*offset="(\d+)"', txt).group(1))
    poff = int(re.search(r'<Points>\s*<DataArray[^>]*offset="(\d+)"', txt, re.S).group(1))
    T = arr(toff, "<f8", npts); P = arr(poff, "<f8", npts*3).reshape(-1, 3)
    return P, T


def topfield(cd):
    P, T = read_vtu(cd); T = T - 273.15
    x, y, z = P[:, 0]*1e3, P[:, 1]*1e3, P[:, 2]*1e3
    m = np.abs(z - 0.745) < 0.15
    return x[m], y[m], T[m]


def case_data(cd):
    P, T = read_vtu(cd); T = T - 273.15
    x, y, z = P[:, 0]*1e3, P[:, 1]*1e3, P[:, 2]*1e3
    out = {}
    for nm, (w, l, cx, cy) in COMPS.items():
        s = (np.abs(x-cx) <= w/2) & (np.abs(y-cy) <= l/2) & (z > 0.7455)
        out[nm] = float(T[s].max()) if s.sum() else np.nan
    out["Tmax"] = float(T.max())
    return out


def fig_matrix():
    vmin, vmax = 25, 130
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 22))
    for i, geo in enumerate(GEOS):
        for j, L in enumerate(LOADS):
            ax = axes[i][j]
            cd = os.path.join(ROOT, "cases", f"{geo}_{LKEY[L]}")
            x, y, T = topfield(cd)
            tp = ax.tricontourf(x, y, T, levels=np.linspace(vmin, vmax, 40), cmap="inferno",
                                vmin=vmin, vmax=vmax, extend="both")
            for nm, (w, l, cx, cy) in COMPS.items():
                ax.add_patch(Rectangle((cx-w/2, cy-l/2), w, l, fill=False, ec="cyan", lw=1.2))
                s = (np.abs(x-cx) <= w/2) & (np.abs(y-cy) <= l/2)
                tm = np.nanmax(T[s]) if s.sum() else np.nan
                ax.text(cx, cy, f"{nm}\n{tm:.0f}", color="white", fontsize=7,
                        ha="center", va="center", weight="bold")
            ax.set_xlim(0, 26.5); ax.set_ylim(0, 50.5); ax.set_aspect("equal")
            ax.set_title(f"{GLAB[geo]} / {L:.1f} A", fontsize=11, weight="bold")
            if i == 2:
                ax.set_xlabel("x (mm)")
            if j == 0:
                ax.set_ylabel("y (mm)")
    fig.suptitle("TPS62933 board top-surface temperature — 3 geometries x 3 loads "
                 "(Tj(U21) analytic on board)", fontsize=14, weight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 0.92, 0.99])
    cax = fig.add_axes([0.94, 0.25, 0.015, 0.5])
    cb = fig.colorbar(tp, cax=cax); cb.set_label("Temperature (degC)")
    fig.savefig(os.path.join(ROOT, "fig_thermal_3x3_matrix.png"), dpi=130)
    print("saved fig_thermal_3x3_matrix.png")


def fig_compare(n):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    xs = [0.6, 1.5, 3.0]
    col = {"fullcu": "tab:blue", "topcut": "tab:orange", "dualcut": "tab:red"}
    for geo in GEOS:
        bm = [n[f"M_{geo}_{LKEY[L]}"]["board_mean_C"] for L in [0.6, 1.5, 3.0]]
        tj = [n[f"M_{geo}_{LKEY[L]}"]["U21_Tj_conservative_C"] for L in [0.6, 1.5, 3.0]]
        l2 = [n[f"M_{geo}_{LKEY[L]}"]["L2_body_max_C"] for L in [0.6, 1.5, 3.0]]
        ax[0].plot(xs, bm, "-o", color=col[geo], label=GLAB[geo])
        ax[1].plot(xs, tj, "-o", color=col[geo], label=GLAB[geo])
        ax[1].plot(xs, l2, "--s", color=col[geo], alpha=0.6)
    ax[0].set_title("Board mean temperature vs load"); ax[1].set_title(
        "U21 junction Tj (solid) & L2 body max (dashed) vs load")
    for a in ax:
        a.set_xlabel("Iout (A)"); a.set_ylabel("Temperature (degC)")
        a.grid(alpha=0.3); a.legend(fontsize=8)
    ax[1].axhline(150, color="red", ls="--", lw=1); ax[1].text(0.62, 151, "150 degC limit", color="red", fontsize=8)
    fig.suptitle("TPS62933 geometry comparison, Tamb=25 degC, natural conv + radiation", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "fig_thermal_3x3_compare.png"), dpi=150)
    print("saved fig_thermal_3x3_compare.png")


def per_case_figs():
    for geo in GEOS:
        for L in LOADS:
            cd = os.path.join(ROOT, "cases", f"{geo}_{LKEY[L]}")
            x, y, T = topfield(cd)
            fig, ax = plt.subplots(figsize=(5.5, 9))
            tp = ax.tricontourf(x, y, T, levels=30, cmap="inferno")
            plt.colorbar(tp, ax=ax).set_label("Temperature (degC)")
            for nm, (w, l, cx, cy) in COMPS.items():
                ax.add_patch(Rectangle((cx-w/2, cy-l/2), w, l, fill=False, ec="cyan", lw=1.3))
                s = (np.abs(x-cx) <= w/2) & (np.abs(y-cy) <= l/2)
                tm = np.nanmax(T[s]) if s.sum() else np.nan
                ax.text(cx, cy, f"{nm}\n{tm:.0f}", color="white", fontsize=8,
                        ha="center", va="center", weight="bold")
            ax.set_xlim(0, 26.5); ax.set_ylim(0, 50.5); ax.set_aspect("equal")
            ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
            ax.set_title(f"{GLAB[geo]} / {L:.1f} A — top-surface T")
            fig.tight_layout()
            fig.savefig(os.path.join(ROOT, f"fig_field_{geo}_{LKEY[L]}.png"), dpi=130)
            plt.close(fig)
    print("saved per-case field figures")


if __name__ == "__main__":
    n = json.load(open(os.path.join(ROOT, "numbers_thermal_3x3.json")))
    fig_matrix()
    fig_compare(n)
    per_case_figs()
    print("POST DONE")

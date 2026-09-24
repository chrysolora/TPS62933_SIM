#!/usr/bin/env python3
"""post_3x3_real.py — figures + comparison for the REAL-polygon rerun.
Outputs (into thermal_3x3_real/):
  fig_thermal_3x3_real_matrix.png     3 geometries x 3 loads top-surface field
  fig_thermal_3x3_rect_vs_real.png    old rectangle vs real polygon comparison
  heat_sources_3x3_real.json          loss budget table
Reads old numbers read-only from ../thermal_3x3/numbers_thermal_3x3.json.
"""
import os, json, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = "/mnt/raid10/sim-work/tps62933/thermal_3x3_real"
OLD = "/mnt/raid10/sim-work/tps62933/thermal_3x3"
os.chdir(ROOT)
COMPS = {"U21": (2.10, 1.60, 8.25, 26.88), "D1": (4.30, 2.60, 14.10, 6.56),
         "L2": (6.50, 6.50, 6.22, 33.74), "L1": (4.20, 4.20, 16.38, 18.75)}
GEOS = ["fullcu", "topcut", "dualcut"]
GLAB = {"fullcu": "FullCopper", "topcut": "TopCutout", "dualcut": "DualCutout"}
LOADS = [3.0, 1.5, 0.6]
LKEY = {3.0: "L30", 1.5: "L15", 0.6: "L06"}
RECT_AREA = {"L2": 6.60 * 7.00, "L1": 4.35 * 4.56}


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
        ln = int(np.frombuffer(raw[p:p + 4], dtype="<u4")[0])
        return np.frombuffer(raw[p + 4:p + 4 + ln], dtype=dt)[:n]
    toff = int(re.search(r'Name="temperature"[^>]*offset="(\d+)"', txt).group(1))
    poff = int(re.search(r'<Points>\s*<DataArray[^>]*offset="(\d+)"', txt, re.S).group(1))
    T = arr(toff, "<f8", npts); P = arr(poff, "<f8", npts * 3).reshape(-1, 3)
    return P, T


def topfield(cd):
    P, T = read_vtu(cd); T = T - 273.15
    x, y, z = P[:, 0] * 1e3, P[:, 1] * 1e3, P[:, 2] * 1e3
    m = np.abs(z - 0.745) < 0.15
    return x[m], y[m], T[m]


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
                ax.add_patch(Rectangle((cx - w / 2, cy - l / 2), w, l, fill=False, ec="cyan", lw=1.2))
                s = (np.abs(x - cx) <= w / 2) & (np.abs(y - cy) <= l / 2)
                tm = np.nanmax(T[s]) if s.sum() else np.nan
                ax.text(cx, cy, f"{nm}\n{tm:.0f}", color="white", fontsize=7,
                        ha="center", va="center", weight="bold")
            # outline real cutout polygons
            for lbl, poly in POLY.get(geo, {}).items():
                px = [p[0] for p in poly] + [poly[0][0]]
                py = [p[1] for p in poly] + [poly[0][1]]
                ax.plot(px, py, "-", color="lime", lw=1.0, alpha=0.9)
            ax.set_xlim(0, 26.5); ax.set_ylim(0, 50.5); ax.set_aspect("equal")
            ax.set_title(f"{GLAB[geo]} / {L:.1f} A", fontsize=11, weight="bold")
            if i == 2:
                ax.set_xlabel("x (mm)")
            if j == 0:
                ax.set_ylabel("y (mm)")
    fig.suptitle("TPS62933 board top-surface temperature - REAL source cutout polygons "
                 "(green = cutout outline)", fontsize=14, weight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 0.92, 0.99])
    cax = fig.add_axes([0.94, 0.25, 0.015, 0.5])
    cb = fig.colorbar(tp, cax=cax); cb.set_label("Temperature (degC)")
    fig.savefig(os.path.join(ROOT, "fig_thermal_3x3_real_matrix.png"), dpi=130)
    print("saved fig_thermal_3x3_real_matrix.png")


def fig_rect_vs_real(nnew, nold):
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    xs = [0.6, 1.5, 3.0]
    col = {"fullcu": "tab:blue", "topcut": "tab:orange", "dualcut": "tab:red"}
    # A: L2 body max
    for geo in GEOS:
        a = [nnew[f"M_{geo}_{LKEY[L]}"]["L2_body_max_C"] for L in xs]
        b = [nold[f"M_{geo}_{LKEY[L]}"]["L2_body_max_C"] for L in xs]
        ax[0][0].plot(xs, a, "-o", color=col[geo], label=f"{GLAB[geo]} real")
        ax[0][0].plot(xs, b, "--s", color=col[geo], alpha=0.55, label=f"{GLAB[geo]} rect(old)")
    ax[0][0].set_title("L2 body max temperature (solid=real polygon, dashed=old rect)")
    ax[0][0].set_ylabel("degC")
    # B: board mean
    for geo in GEOS:
        a = [nnew[f"M_{geo}_{LKEY[L]}"]["board_mean_C"] for L in xs]
        b = [nold[f"M_{geo}_{LKEY[L]}"]["board_mean_C"] for L in xs]
        ax[0][1].plot(xs, a, "-o", color=col[geo], label=f"{GLAB[geo]} real")
        ax[0][1].plot(xs, b, "--s", color=col[geo], alpha=0.55)
    ax[0][1].set_title("Board mean temperature")
    ax[0][1].set_ylabel("degC")
    # C: cutout area real vs rect
    labels = ["cut L2\n(top)", "cut L1\n(top)", "cut L2\n(dual)", "cut L1\n(dual)"]
    real = [POLY["topcut"]["L2"], POLY["topcut"]["L1"], POLY["dualcut"]["L2"], POLY["dualcut"]["L1"]]
    realA = []
    for p in real:
        a = 0.0
        for i in range(len(p)):
            x1, y1 = p[i]; x2, y2 = p[(i + 1) % len(p)]
            a += x1 * y2 - x2 * y1
        realA.append(abs(a) / 2.0)
    rectA = [RECT_AREA["L2"], RECT_AREA["L1"], RECT_AREA["L2"], RECT_AREA["L1"]]
    xx = np.arange(4); ww = 0.38
    ax[1][0].bar(xx - ww / 2, rectA, ww, label="old rect (bbox)", color="tab:gray")
    ax[1][0].bar(xx + ww / 2, realA, ww, label="real polygon", color="tab:green")
    for i in range(4):
        ax[1][0].text(i, max(rectA[i], realA[i]) + 0.6,
                      f"-{100*(rectA[i]-realA[i])/rectA[i]:.0f}%", ha="center", fontsize=8)
    ax[1][0].set_xticks(xx); ax[1][0].set_xticklabels(labels, fontsize=8)
    ax[1][0].set_ylabel("area (mm^2)"); ax[1][0].legend(fontsize=8)
    ax[1][0].set_title("Cutout area: old rectangle vs real polygon (over-cut %)")
    # D: cutout effect on L2 at 3A  (cut - fullcu)
    geos2 = ["topcut", "dualcut"]
    dreal = [nnew[f"M_{g}_L30"]["L2_body_max_C"] - nnew["M_fullcu_L30"]["L2_body_max_C"] for g in geos2]
    dold = [nold[f"M_{g}_L30"]["L2_body_max_C"] - nold["M_fullcu_L30"]["L2_body_max_C"] for g in geos2]
    xx = np.arange(2); ww = 0.38
    ax[1][1].bar(xx - ww / 2, dold, ww, label="old rect", color="tab:gray")
    ax[1][1].bar(xx + ww / 2, dreal, ww, label="real polygon", color="tab:green")
    for i in range(2):
        ax[1][1].text(i - ww / 2, dold[i] + 0.1, f"{dold[i]:.2f}", ha="center", fontsize=8)
        ax[1][1].text(i + ww / 2, dreal[i] + 0.1, f"{dreal[i]:.2f}", ha="center", fontsize=8)
    ax[1][1].set_xticks(xx); ax[1][1].set_xticklabels([GLAB[g] for g in geos2])
    ax[1][1].set_ylabel("dT (degC)")
    ax[1][1].legend(fontsize=8)
    ax[1][1].set_title("Cutout-induced L2 hotspot at 3.0 A (L2 - FullCopper)")
    for a in (ax[1][1],):
        a.grid(alpha=0.3)
    for a in (ax[0][0], ax[0][1]):
        a.set_xlabel("Iout (A)"); a.grid(alpha=0.3)
        a.legend(fontsize=7, ncol=2)
    fig.suptitle("TPS62933 thermal 3x3: OLD assumed rectangle vs REAL source cutout polygon",
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(ROOT, "fig_thermal_3x3_rect_vs_real.png"), dpi=140)
    print("saved fig_thermal_3x3_rect_vs_real.png")


def heat_sources():
    import importlib.util
    spec = importlib.util.spec_from_file_location("mc", os.path.join(ROOT, "make_case_3x3.py"))
    mc = importlib.util.module_from_spec(spec); spec.loader.exec_module(mc)
    out = {"note": "Loss budget identical to thermal_final/thermal_3x3 (no residual term). "
                   "Only the keep-out GEOMETRY changed to real source polygons.",
           "loads": {}}
    for Iout in [3.0, 1.5, 0.6]:
        p, Iin = mc.losses(Iout)
        P = {"D1": p["D1"], "L2": p["L2"], "L1": p["L1"],
             "U21_cond": p["U21_cond"], "U21_sw": p["U21_sw"]}
        P["U21_total"] = P["U21_cond"] + P["U21_sw"]
        P["total"] = P["D1"] + P["L2"] + P["L1"] + P["U21_total"]
        P["Iin_A"] = Iin; P["Iout_A"] = Iout
        out["loads"][f"{Iout}A"] = {k: round(v, 5) for k, v in P.items()}
    json.dump(out, open(os.path.join(ROOT, "heat_sources_3x3_real.json"), "w"), indent=1)
    print("saved heat_sources_3x3_real.json")


POLY = {}
def load_poly():
    global POLY
    d = json.load(open(os.path.join(ROOT, "cutout_polygons.json")))
    for g in ("topcut", "dualcut"):
        POLY[g] = {lbl: rec["vertices_mesh_mm"] for lbl, rec in d["polygons"][g].items()}


if __name__ == "__main__":
    load_poly()
    nnew = json.load(open(os.path.join(ROOT, "numbers_thermal_3x3_real.json")))
    nold = json.load(open(os.path.join(OLD, "numbers_thermal_3x3.json")))
    fig_matrix()
    fig_rect_vs_real(nnew, nold)
    heat_sources()
    print("POST DONE")
